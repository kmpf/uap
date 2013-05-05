'''
This module can be used to launch child processes and wait for them.
Processes may either run on their own or pipelines can be built with them.
'''

import sys
sys.path.append('./include/steps')
import copy
import datetime
import fcntl
import hashlib
import os
import signal
import subprocess
import tempfile
import time
import yaml


class UnixPipeline(object):

    TAIL_LENGTH = 1024
    COPY_BLOCK_SIZE = 4194304
    SIGTERM_TIMEOUT = 5

    # signal names for numbers... kudos to http://stackoverflow.com/questions/2549939/get-signal-names-from-numbers-in-python
    SIGNAL_NAMES = dict((getattr(signal, n), n) for n in dir(signal) if n.startswith('SIG') and '_' not in n )

    class Pipeline(object):
        '''
        This class can be used to chain multiple processes together.
        '''
        def __init__(self):
            self.append_calls = []

        def append(self, args, stdout_path = None, stderr_path = None):
            call = {
                'args': copy.deepcopy(args),
                'stdout_path': copy.copy(stdout_path),
                'stderr_path': copy.copy(stderr_path)
            }
            self.append_calls.append(call)
            
    class TimeoutException(Exception):
        pass

    def timeout_handler(signum, frame):
        raise TimeoutException()
    
    def __init__(self, step):
        # the current step this UnixPipeline is used in (for temporary paths etc.)
        self.step = step
        
        # log entries
        self.log_entries = []
        
        # dict of PID -> path to copy process report
        self.copy_process_reports = {}
        
        # set of currently running PIDs
        # whenever a child process exits, its PID will get removed from this
        self.running_procs = set()
        
        # list of PIDs in the order that processes were launched
        # whenever a child process exits, its PID will remain in here
        self.proc_order = []
        
        # dict of PID -> process info
        # whenever a child process exits, its PID will remain in here
        self.proc_details = {}
        
        # list of temp paths to clean up
        self.temp_paths = []
        
        # list of commands to be launched
        self.launch_calls = []
        
    def __enter__(self):
        return self
        
    def __exit__(self, type, value, traceback):
        # now launch all processes...
        self._launch_all_processes()

        # ...and wait until all child processes have exited
        try:
            self._wait()
        except:
            # pass log to step even if there was a problem
            self.step.append_pipeline_log(self.get_log())
            raise
        
        # if there was no exception, still pass log to step
        self.step.append_pipeline_log(self.get_log())
        
        # remove all temporary files we know of
        for _ in self.temp_paths:
            try:
                os.unlink(_)
            except OSError:
                pass
        
    def get_temporary_path(self, prefix, designation = None):
        path = self.step.get_temporary_path(prefix, designation)
        self.temp_paths.append(path)
        return path
        
    def get_temporary_fifo(self, prefix, designation = None):
        path = self.step.get_temporary_fifo(prefix, designation)
        self.temp_paths.append(path)
        return path
    
    def create_pipeline(self):
        p = self.Pipeline()
        self.launch_calls.append(p)
        return p
    
    def launch(self, args, stdout_path = None, stderr_path = None):
        call = {
            'args': copy.deepcopy(args),
            'stdout_path': copy.copy(stdout_path),
            'stderr_path': copy.copy(stderr_path)
        }
        
        self.launch_calls.append(call)
        
    def log(self, message):
        '''
        Append a message to the pipeline log.
        '''
        formatted_message = "[%s] %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message)
        self.log_entries.append(formatted_message)
        
    def get_log(self):
        '''
        Return the log as a dictionary.
        '''
        log = dict()
        log['processes'] = []
        proc_details_copy = dict()
        for pid in self.proc_order:
            if 'listener_for' in self.proc_details[pid]:
                attach_info = self.proc_details[pid]['listener_for']
                proc_details_copy[attach_info[0]][attach_info[1] + '_copy'] = copy.deepcopy(self.proc_details[pid])
                del proc_details_copy[attach_info[0]][attach_info[1] + '_copy']['listener_for']
            else:
                proc_details_copy[pid] = copy.deepcopy(self.proc_details[pid])
            
        for pid in self.proc_order:
            if pid in proc_details_copy:
                log['processes'].append(proc_details_copy[pid])
                
        log['log'] = copy.deepcopy(self.log_entries)
        
        return log
        
    def _launch_all_processes(self):
        for info in self.launch_calls:
            if info.__class__ == UnixPipeline.Pipeline:
                pipeline = info
                use_stdin = None
                last_pid = None
                for index, info in enumerate(pipeline.append_calls):
                    use_stdin, pid = self._do_launch(info, 
                        index < len(pipeline.append_calls) - 1, use_stdin)
                    if last_pid is not None:
                        self.proc_details[pid]['use_stdin_of'] = last_pid
                    last_pid = pid
            else:
                self._do_launch(info)
    
    def _do_launch(self, info, keep_stdout_open = False, use_stdin = None):
        args = copy.deepcopy(info['args'])
        stdout_path = copy.copy(info['stdout_path'])
        stderr_path = copy.copy(info['stderr_path'])
        
        program_name = copy.deepcopy(args[0])
        if program_name.__class__ == list:
            new_args = args[0]
            program_name = new_args[-1]
            new_args.extend(args[1:])
            args = new_args

        # launch the process and always pipe stdout and stderr because we
        # want to watch both streams, regardless of whether stdout should 
        # be passed on to another process
        proc = subprocess.Popen(args,
            stdin = use_stdin,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            preexec_fn = os.setsid
        )
        pid = proc.pid
        print("%d launched: %s" % (pid, args))

        self.running_procs.add(pid)
        self.proc_order.append(pid)
        self.proc_details[pid] = {
            'name': os.path.basename(program_name),
            'start_time': datetime.datetime.now(),
            'args': args,
            'pid': pid
        }
        self.log("Launched %s as PID %d." % (' '.join(args), pid))

        pipe = None
        if keep_stdout_open:
            pipe = os.pipe()
        
        for which in ['stdout', 'stderr']:
            report_path = self.get_temporary_path("%s-report" % which)
            sink_path = stdout_path if which == 'stdout' else stderr_path
            listener_pid = self._do_launch_copy_process(
                proc.stdout if which == 'stdout' else proc.stderr,
                sink_path,
                report_path, pid, which, 
                pipe if which == 'stdout' else None)
            
            self.copy_process_reports[listener_pid] = report_path
            
            if sink_path is not None:
                self.proc_details[listener_pid]['sink'] = os.path.basename(sink_path)
                self.proc_details[listener_pid]['sink_full_path'] = sink_path

        if keep_stdout_open:
            os.close(pipe[1])
            return pipe[0], pid
        else:
            return None, pid

    def _do_launch_copy_process(self, fin, fout_path, report_path, parent_pid, which, pipe):
        pid = os.fork()
        if pid == 0:
            if pipe is not None:
                os.close(pipe[0])
            fdout = None
            if fout_path is not None:
                fdout = os.open(fout_path, os.O_WRONLY|os.O_CREAT|os.O_TRUNC)
                
            checksum = hashlib.sha1()
            tail = ''
            length = 0
            newline_count = 0
            
            while True:
                block = fin.read(UnixPipeline.COPY_BLOCK_SIZE)
                if len(block) == 0:
                    # fin has been closed, let's call it a day
                    break
                    
                # update checksum
                checksum.update(block)
                
                # update tail
                if len(block) >= UnixPipeline.TAIL_LENGTH:
                    tail = block[-UnixPipeline.TAIL_LENGTH:]
                else:
                    keep_length = UnixPipeline.TAIL_LENGTH - len(block)
                    tail = tail[0:keep_length] + block
                    
                # update length
                length += len(block)
                
                # update newline_count
                newline_count += block.count('\n')
                
                # write block to output file
                if fdout is not None:
                    bytes_written = os.write(fdout, block)
                    if bytes_written != len(block):
                        sys.stderr.write("Could not write to fdout.\n")
                        os._exit(1)
                    
                # write block to pipe
                if pipe is not None:
                    bytes_written = os.write(pipe[1], block)
                    if bytes_written != len(block):
                        sys.stderr.write("Could not write to pipe.\n")
                        os._exit(2)
                    
            # we're finished, close everything
            if fdout is not None:
                os.close(fdout)
            if pipe is not None:
                os.close(pipe[1])

            # write report
            with open(report_path, 'w') as freport:
                report = dict()
                report['sha1'] = checksum.hexdigest()
                report['tail'] = tail
                report['length'] = length
                report['lines'] = newline_count
                freport.write(yaml.dump(report))
                
            os._exit(0)
        else:
            print("%d launched: %s of %d" % (pid, which, parent_pid))
            self.running_procs.add(pid)
            self.proc_order.append(pid)
            self.proc_details[pid] = {
                'name': '[stream listener for %s of PID %d]' % (which, parent_pid),
                'start_time': datetime.datetime.now(),
                'pid': pid,
                'listener_for': [parent_pid, which]
            }
            self.log("Launched a copy process with PID %d to capture %s of PID %d." % (pid, which, parent_pid))
            if fout_path is not None:
                self.log("...which gets also redirected to %s" % fout_path)
            return pid
                
    def _wait(self):
        something_went_wrong = False
        while True:
            try:
                # wait for the next child process to exit
                pid, exit_code_with_signal = os.wait()
                print("%d exited: %d" % (pid, exit_code_with_signal))
                
                # remove pid from self.running_procs
                self.running_procs.remove(pid)
                
                self.proc_details[pid]['end_time'] = datetime.datetime.now()
                
                signal_number = exit_code_with_signal & 255
                exit_code = exit_code_with_signal >> 8
                what_happened = "has exited with exit code %d" % exit_code
                if signal_number > 0:
                    what_happened = "has received signal %d" % signal_number
                    if signal_number in UnixPipeline.SIGNAL_NAMES:
                        what_happened = "has received %s (signal number %d)" % (UnixPipeline.SIGNAL_NAMES[signal_number], signal_number)
                        
                self.log("%s (PID %d) %s." % (self.proc_details[pid]['name'], pid, what_happened))

                if signal_number == 0:
                    self.proc_details[pid]['exit_code'] = exit_code
                else:
                    self.proc_details[pid]['signal'] = signal_number
                    if signal_number in UnixPipeline.SIGNAL_NAMES:
                        self.proc_details[pid]['signal_name'] = UnixPipeline.SIGNAL_NAMES[signal_number]
                        
                try:
                    if pid in self.copy_process_reports:
                        report_path = self.copy_process_reports[pid]
                        if os.path.exists(report_path):
                            report = yaml.load(open(report_path, 'r'))
                            os.unlink(report_path)
                            self.proc_details[pid].update(report)
                except Exception:
                    # swallow possible exceptions because getting the checksum and
                    # the tail is nice to have, but no cause for the entire process
                    # to fail
                    raise
                    pass
                
            except UnixPipeline.TimeoutException:
                log("Timeout, killing all child processes now.")
                self._kill_all_child_processes()
            except OSError:
                # no more children running, we are done
                signal.alarm(0)
                self.log("Cancelling timeout (if there was one), all child processes have exited.")
                break
            else:
                if exit_code_with_signal != 0:
                    # Oops, something went wrong. See what happens and terminate
                    # all child processes in a few seconds.
                    something_went_wrong = True
                    signal.signal(signal.SIGALRM, UnixPipeline.timeout_handler) 
                    self.log("Terminating all children in %d seconds..." % UnixPipeline.SIGTERM_TIMEOUT)
                    signal.alarm(UnixPipeline.SIGTERM_TIMEOUT)
                            
        if something_went_wrong:
            self.log("Pipeline crashed.")
            raise StandardError("Pipeline crashed.")        
        
    def _kill_all_child_processes(self):
        '''
        Kill all child processes launched via this module by sending a SIGTERM to each of them.
        '''
        for pid in self.running_procs:
            try:
                os.kill(pid, signal.SIGTERM)
                log('Sending SIGTERM to process %s (PID %d).' % (self.proc_details[pid]['name'], pid))
            except OSError:
                pass
