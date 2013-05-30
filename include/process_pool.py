'''
This module can be used to launch child processes and wait for them.
Processes may either run on their own or pipelines can be built with them.
'''

import sys
sys.path.append('./include/steps')
import copy
import datetime
import errno
import fcntl
import hashlib
import misc
import os
import psutil
import signal
import subprocess
import tempfile
import time
import traceback
import yaml

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException()

def restore_sigpipe_handler():
    # http://www.chiark.greenend.org.uk/ucgi/~cjwatson/blosxom/2009-07-02-python-sigpipe.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.setsid()
    
class ProcessPool(object):

    TAIL_LENGTH = 1024
    COPY_BLOCK_SIZE = 4194304
    SIGTERM_TIMEOUT = 10
    
    process_watcher_pid = None

    current_instance = None
    process_pool_is_dead = False

    # signal names for numbers... kudos to http://stackoverflow.com/questions/2549939/get-signal-names-from-numbers-in-python
    SIGNAL_NAMES = dict((getattr(signal, n), n) for n in dir(signal) if n.startswith('SIG') and '_' not in n )

    class Pipeline(object):
        '''
        This class can be used to chain multiple processes together.
        '''
        def __init__(self, pool):
            pool.launch_calls.append(self)
            self.append_calls = []
            
        def __enter__(self):
            return self
            
        def __exit__(self, type, value, traceback):
            pass

        def append(self, args, stdout_path = None, stderr_path = None, hints = {}):
            call = {
                'args': copy.deepcopy(args),
                'stdout_path': copy.copy(stdout_path),
                'stderr_path': copy.copy(stderr_path),
                'hints': copy.deepcopy(hints)
            }
            self.append_calls.append(call)
            
    def __init__(self, step):
        if ProcessPool.process_pool_is_dead:
            raise StandardError("We have encountered an error, stopping now...")
        
        # the current step this ProcessPool is used in (for temporary paths etc.)
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

        # we must keep Popen objects around, or their destructor will be called
        self.popen_procs = {}
        
        # dict of PID -> process info
        # whenever a child process exits, its PID will remain in here
        self.proc_details = {}
        
        self.process_watcher_report = dict()
        
        # list of temp paths to clean up
        self.temp_paths = []
        
        # list of commands to be launched
        self.launch_calls = []
        
        # List of processes we killed deliberately. Look: every time a 
        # within a pipeline exits, we SIGTERM its predecessor. This is
        # necessary because otherwise, stuff is hanging forever.
        self.ok_to_fail = set()
        
        self.copy_processes_for_pid = dict()
        
    def __enter__(self):
        if ProcessPool.current_instance is not None:
            raise StandardError("Sorry, only one instance of ProcessPool allowed at a time.")
        ProcessPool.current_instance = self
        
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
            
        ProcessPool.current_instance = None
        
    def get_temporary_path(self, prefix, designation = None):
        path = self.step.get_temporary_path(prefix, designation)
        self.temp_paths.append(path)
        return path
        
    def get_temporary_fifo(self, prefix, designation = None):
        path = self.step.get_temporary_fifo(prefix, designation)
        self.temp_paths.append(path)
        return path
    
    def launch(self, args, stdout_path = None, stderr_path = None, hints = {}):
        call = {
            'args': copy.deepcopy(args),
            'stdout_path': copy.copy(stdout_path),
            'stderr_path': copy.copy(stderr_path),
            'hints': copy.deepcopy(hints)
        }
        
        self.launch_calls.append(call)
        
    def log(self, message):
        '''
        Append a message to the pipeline log.
        '''
        formatted_message = "[%s] %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message)
        #sys.stderr.write(formatted_message + "\n")
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
        log['process_watcher'] = copy.deepcopy(self.process_watcher_report)
        log['ok_to_fail'] = copy.deepcopy(self.ok_to_fail)
        
        return log
        
    def _launch_all_processes(self):
        for info in self.launch_calls:
            if info.__class__ == ProcessPool.Pipeline:
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
        hints = copy.deepcopy(info['hints'])
        
        program_name = copy.deepcopy(args[0])
        if program_name.__class__ == list:
            new_args = args[0]
            program_name = new_args[-1]
            new_args.extend(args[1:])
            args = new_args

        # launch the process and always pipe stdout and stderr because we
        # want to watch both streams, regardless of whether stdout should 
        # be passed on to another process
        proc = subprocess.Popen(
            args,
            stdin = use_stdin,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            preexec_fn = restore_sigpipe_handler,
            close_fds = True
        )
        pid = proc.pid
        self.popen_procs[pid] = proc

        self.running_procs.add(pid)
        self.proc_order.append(pid)
        self.proc_details[pid] = {
            'name': os.path.basename(program_name),
            'start_time': datetime.datetime.now(),
            'args': args,
            'pid': pid,
            'hints': hints
        }
        self.log("Launched %s as PID %d." % (' '.join(args), pid))

        pipe = None
        if keep_stdout_open:
            pipe = os.pipe()
        
        self.copy_processes_for_pid[pid] = list()
        
        for which in ['stdout', 'stderr']:
            report_path = self.get_temporary_path("%s-report" % which)
            sink_path = stdout_path if which == 'stdout' else stderr_path
            listener_pid = self._do_launch_copy_process(
                proc.stdout if which == 'stdout' else proc.stderr,
                sink_path,
                report_path, pid, which, 
                pipe if which == 'stdout' else None)
                
            self.copy_processes_for_pid[pid].append(listener_pid)
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
            
            def write_report_and_exit():
                try:
                    # write report
                    with open(report_path, 'w') as freport:
                        #sys.stderr.write("[%d] writing report...\n" % os.getpid())
                        report = dict()
                        report['sha1'] = checksum.hexdigest()
                        report['tail'] = tail
                        report['length'] = length
                        report['lines'] = newline_count
                        freport.write(yaml.dump(report))
                except:
                    pass
                
                os._exit(0)
                
            def sigpipe_handler(signum, frame):
                write_report_and_exit()
            
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
            signal.signal(signal.SIGINT, signal.SIG_DFL)
            signal.signal(signal.SIGPIPE, sigpipe_handler)
            os.setsid()
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
                #sys.stderr.write("[%d] reading from fin\n" % os.getpid())
                block = fin.read(ProcessPool.COPY_BLOCK_SIZE)
                #sys.stderr.write("[%d] actually read %d bytes from fin\n" % (os.getpid(), len(block)))
                if len(block) == 0:
                    # fin reports EOF, let's call it a day
                    #sys.stderr.write("[%d] fin is at EOF\n" % os.getpid())
                    break
                    
                # update checksum
                checksum.update(block)
                
                # update tail
                if len(block) >= ProcessPool.TAIL_LENGTH:
                    tail = block[-ProcessPool.TAIL_LENGTH:]
                else:
                    keep_length = ProcessPool.TAIL_LENGTH - len(block)
                    tail = tail[0:keep_length] + block
                    
                # update length
                length += len(block)
                
                # update newline_count
                newline_count += block.count('\n')
                
                # write block to output file
                if fdout is not None:
                    #sys.stderr.write("[%d] writing %d bytes to fdout\n" % (os.getpid(), len(block)))
                    bytes_written = os.write(fdout, block)
                    #sys.stderr.write("[%d] actually wrote %d bytes to fdout\n" % (os.getpid(), bytes_written))
                    if bytes_written != len(block):
                        #sys.stderr.write("Could not write to fdout.\n")
                        os._exit(1)
                    
                # write block to pipe
                if pipe is not None:
                    #sys.stderr.write("[%d] writing %d bytes to pipe[1]\n" % (os.getpid(), len(block)))
                    bytes_written = os.write(pipe[1], block)
                    #sys.stderr.write("[%d] actually wrote %d bytes to pipe[1]\n" % (os.getpid(), bytes_written))
                    if bytes_written != len(block):
                        #sys.stderr.write("Could not write to pipe.\n")
                        os._exit(2)
                    
            # we're finished, close everything
            #sys.stderr.write("[%d] closing fin...\n" % os.getpid())
            fin.close()
            #sys.stderr.write("[%d] done closing fin...\n" % os.getpid())
            if fdout is not None:
                #sys.stderr.write("[%d] closing fdout...\n" % os.getpid())
                os.close(fdout)
                #sys.stderr.write("[%d] done closing fdout...\n" % os.getpid())
            if pipe is not None:
                #sys.stderr.write("[%d] closing pipe[1]...\n" % os.getpid())
                os.close(pipe[1])
                #sys.stderr.write("[%d] done closing pipe[1]...\n" % os.getpid())

            write_report_and_exit()
                
            os._exit(0)
        else:
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
        self.log("Now launching process watcher and waiting for all child processes to exit.")
        watcher_report_path = self.step.get_temporary_path('watcher-report')
        watcher_pid = self._launch_process_watcher(watcher_report_path)
        ProcessPool.process_watcher_pid = watcher_pid
        something_went_wrong = False
        while True:
            if len(self.running_procs) == 0:
                break
            try:
                # wait for the next child process to exit
                pid, exit_code_with_signal = os.wait()
                if pid == watcher_pid:
                    ProcessPool.process_watcher_pid = None
                    try:
                        self.process_watcher_report = yaml.load(open(watcher_report_path))
                        os.unlink(watcher_report_path)
                    except:
                        print("Warning: Couldn't load watcher report from %s." % watcher_report_path)
                        pass
                    # the process watcher has terminated, which is cool, I guess
                    # (if it's the last child process, anyway)
                    continue

                try:
                    # remove pid from self.running_procs
                    self.running_procs.remove(pid)
                except:
                    if pid != os.getpid():
                        #raise StandardError("Caught a process which we didn't know: %d." % pid)
                        sys.stderr.write("Note: Caught a process which we didn't know: %d.\n" % pid)

                if pid in self.proc_details:
                    self.proc_details[pid]['end_time'] = datetime.datetime.now()
                
                signal_number = exit_code_with_signal & 255
                exit_code = exit_code_with_signal >> 8
                what_happened = "has exited with exit code %d" % exit_code
                if signal_number > 0:
                    what_happened = "has received signal %d" % signal_number
                    if signal_number in ProcessPool.SIGNAL_NAMES:
                        what_happened = "has received %s (signal number %d)" % (ProcessPool.SIGNAL_NAMES[signal_number], signal_number)
                        
                if pid in self.proc_details:
                    self.log("%s (PID %d) %s." % (self.proc_details[pid]['name'], pid, what_happened))
                else:
                    self.log("PID %d %s." % (pid, what_happened))

                if pid in self.proc_details:
                    if signal_number == 0:
                        self.proc_details[pid]['exit_code'] = exit_code
                    else:
                        self.proc_details[pid]['signal'] = signal_number
                        if signal_number in ProcessPool.SIGNAL_NAMES:
                            self.proc_details[pid]['signal_name'] = ProcessPool.SIGNAL_NAMES[signal_number]
                            
                    # now kill it's preceding process, if this is from a pipeline
                    if 'use_stdin_of' in self.proc_details[pid]:
                        pidlist = list()
                        pidlist.append(self.copy_processes_for_pid[self.proc_details[pid]['use_stdin_of']][0])
                        pidlist.append(self.copy_processes_for_pid[self.proc_details[pid]['use_stdin_of']][1])
                        pidlist.append(self.proc_details[pid]['use_stdin_of'])
                        for kpid in pidlist:
                            self.log("Now killing %d, the predecessor of %d." % (kpid, pid))
                            self.ok_to_fail.add(kpid)
                            try:
                                os.kill(kpid, signal.SIGPIPE)
                            except OSError, e:
                                if e.errno == errno.ESRCH:
                                    self.log("Couldn't kill %d: no such process." % kpid)
                                    pass
                                else:
                                    raise
                                
                    if pid in self.copy_process_reports:
                        report_path = self.copy_process_reports[pid]
                        if os.path.exists(report_path):
                            report = yaml.load(open(report_path, 'r'))
                            os.unlink(report_path)
                            if report is not None:
                                self.proc_details[pid].update(report)
                
            except TimeoutException:
                self.log("Timeout, killing all child processes now.")
                ProcessPool.kill_all_child_processes()
            except OSError, e:
                if e.errno == errno.ECHILD:
                    # no more children running, we are done
                    sys.stderr.write("ProcessPool: There are no child processes left, exiting.\n")
                    signal.alarm(0)
                    self.log("Cancelling timeout (if there was one), all child processes have exited.")
                    break
                elif e.errno == errno.EINTR:
                    # a system call was interrupted, pfft.
                    pass
                else:
                    raise
            else:
                if exit_code_with_signal != 0:
                    if not pid in self.ok_to_fail:
                        # Oops, something went wrong. See what happens and terminate
                        # all child processes in a few seconds.
                        something_went_wrong = True
                        signal.signal(signal.SIGALRM, timeout_handler)
                        self.log("Terminating all children in %d seconds..." % ProcessPool.SIGTERM_TIMEOUT)
                        signal.alarm(ProcessPool.SIGTERM_TIMEOUT)
                        
        # now wait for the watcher process, if it still exists
        try:
            os.waitpid(watcher_pid, 0)
            try:
                self.process_watcher_report = yaml.load(open(watcher_report_path))
                os.unlink(watcher_report_path)
            except:
                print("Warning: Couldn't load watcher report from %s." % watcher_report_path)
                pass
        except OSError, e:
            if e.errno == errno.ESRCH:
                pass
            elif e.errno == errno.ECHILD:
                pass
            else:
                raise
                            
        if something_went_wrong:
            self.log("Pipeline crashed.")
            raise StandardError("Pipeline crashed.")        
        
    def _launch_process_watcher(self, watcher_report_path):
        super_pid = os.getpid()
        
        watcher_pid = os.fork()
        if watcher_pid == 0:
            try:
                signal.signal(signal.SIGTERM, signal.SIG_DFL)
                signal.signal(signal.SIGINT, signal.SIG_IGN)
                called_cpu_stat_for_childpid = set()
                procs = {}
                procs[super_pid] = psutil.Process(super_pid)
                procs[os.getpid()] = psutil.Process(os.getpid())
                for pid in self.running_procs:
                    try:
                        procs[pid] = psutil.Process(pid)
                    except psutil._error.NoSuchProcess:
                        pass

                pid_list = copy.deepcopy(procs.keys())
                for pid in pid_list:
                    proc = procs[pid]
                    try:
                        cpu_percent = proc.get_cpu_percent(interval = None)
                        # add values for all children
                        if pid != super_pid:
                            for p in proc.get_children(recursive = True):
                                try:
                                    cpu_percent = p.get_cpu_percent(interval = None)
                                    called_cpu_stat_for_childpid.add(p.pid)
                                except psutil._error.NoSuchProcess:
                                    pass
                    except psutil._error.NoSuchProcess:
                        del procs[pid]

                time.sleep(0.1)

                iterations = 0
                delay = 0.5
                max_data = dict()
                while True:
                    pid_list = copy.deepcopy(procs.keys())
                    sum_data = dict()
                    for pid in pid_list:
                        proc = procs[pid]
                        try:
                            data = dict()
                            data['cpu_percent'] = proc.get_cpu_percent(interval = None)
                            data['memory_percent'] = proc.get_memory_percent()
                            memory_info = proc.get_memory_info()
                            data['rss'] = memory_info.rss
                            data['vms'] = memory_info.vms
                            
                            # add values for all children
                            if pid != super_pid:
                                for p in proc.get_children(recursive = True):
                                    try:
                                        v = p.get_cpu_percent(interval = None)
                                        if p.pid in called_cpu_stat_for_childpid:
                                            data['cpu_percent'] += v
                                        called_cpu_stat_for_childpid.add(p.pid)
                                        data['memory_percent'] += p.get_memory_percent()
                                        memory_info = p.get_memory_info()
                                        data['rss'] += memory_info.rss
                                        data['vms'] += memory_info.vms
                                    except psutil._error.NoSuchProcess:
                                        pass
                            
                            if not pid in max_data:
                                max_data[pid] = copy.deepcopy(data)
                            for k, v in data.items():
                                max_data[pid][k] = max(max_data[pid][k], v)

                            if len(sum_data) == 0:
                                for k, v in data.items():
                                    sum_data[k] = 0.0

                            for k, v in data.items():
                                sum_data[k] += v

                        except psutil._error.NoSuchProcess:
                            del procs[pid]

                    if not 'sum' in max_data:
                        max_data['sum'] = copy.deepcopy(sum_data)
                    for k, v in sum_data.items():
                        max_data['sum'][k] = max(max_data['sum'][k], v)

                    if len(procs) <= 2:
                        # there's nothing more to watch, write report and exit
                        # (now there's only the controlling python process and
                        # the process watcher itself
                        with open(watcher_report_path, 'w') as f:
                            report = dict()
                            report['max'] = max_data
                            f.write(yaml.dump(report, default_flow_style = False))

                        os._exit(0)

                    iterations += 1
                    if iterations == 10:
                        delay = 1
                    if iterations == 20:
                        delay = 2
                    if iterations == 30:
                        delay = 10
                    time.sleep(delay)
            finally:
                os._exit(0)
        else:
            return watcher_pid

    @classmethod
    def kill(cls):
        '''
        Kills all user-launched processes. After that, the remaining process will end 
        and a report will be written.
        '''
        ProcessPool.process_pool_is_dead = True
        
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(ProcessPool.SIGTERM_TIMEOUT)
        if ProcessPool.current_instance is not None:
            for pid in ProcessPool.current_instance.copy_processes_for_pid.keys():
                try:
                    os.kill(pid, signal.SIGTERM)
                except:
                    pass

    @classmethod
    def kill_all_child_processes(cls):
        '''
        Kill all child processes of this process by sending a SIGTERM to each of them.
        This includes all children which were not launched by this module, and
        their children etc.
        '''
        proc = psutil.Process(os.getpid())
        for p in proc.get_children(recursive = True):
            try:
                p.terminate()
            except psutil._error.NoSuchProcess:
                pass
        