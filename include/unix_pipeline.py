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

TAIL_LENGTH = 1024
COPY_BLOCK_SIZE = 4194304
SIGTERM_TIMEOUT = 5

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException()

# for better error messages, this dict holds the process names for PIDs
name_for_pid = {}

sha1_checksum_for_file_basename = {}

up_log = []
copy_process_reports = {}
seal_these_pipeline_instances = []

proc_order = []
proc_details = {}

def clear():
    name_for_pid = dict()
    sha1_checksum_for_file_basename = dict()
    up_log = list()
    copy_process_reports = dict()
    seal_these_pipeline_instances = list()
    proc_order = list()
    proc_details = dict()
    
def add_proc_info(pid, name = None, args = None):
    proc_order.append(pid)
    proc_details[pid] = {}
    proc_details[pid]['pid'] = pid
    if name is not None:
        proc_details[pid]['name'] = name
    if args is not None:
        proc_details[pid]['args'] = ' '.join(args)

def log(message):
    '''
    Append a message to the pipeline log.
    '''
    #sys.stderr.write(message + "\n")
    up_log.append(message)
    
def get_log():
    '''
    Return the log as a dictionary.
    '''
    log = dict()
    log['processes'] = []
    proc_details_copy = dict()
    for pid in proc_order:
        if 'attach_to' in proc_details[pid]:
            attach_info = proc_details[pid]['attach_to']
            proc_details_copy[attach_info[0]][attach_info[1] + '_copy'] = copy.deepcopy(proc_details[pid])
            del proc_details_copy[attach_info[0]][attach_info[1] + '_copy']['attach_to']
        else:
            proc_details_copy[pid] = copy.deepcopy(proc_details[pid])
        
    for pid in proc_order:
        if pid in proc_details_copy:
            log['processes'].append(proc_details_copy[pid])
            
    log['log'] = up_log
    return log
    
def mkfifo(suffix = ''):
    '''
    Create a temporary FIFO and return its path.
    '''
    _, path = tempfile.mkstemp(suffix)
    os.close(_)
    os.unlink(path)
    os.mkfifo(path)
    return path

def temppath(suffix = ''):
    '''
    Return a temporary file path.
    '''
    _, path = tempfile.mkstemp(suffix)
    os.close(_)
    os.unlink(path)
    return path

def kill_all_child_processes():
    '''
    Kill all child processes launched via this module by sending a SIGTERM to each of them.
    '''
    for pid, name in name_for_pid.items():
        try:
            os.kill(pid, signal.SIGTERM)
            log('Sending SIGTERM to process %s (PID %d).' % (name, pid))
            sys.stderr.write("Sending SIGTERM to " + name + " (PID " + str(pid) + ").\n")
            sys.stderr.flush()
        except OSError:
            pass
        

def launch_copy_process(fin, fout, report_path, other_pid, which, pipe, use_pipe):
    '''
    Launch a copy process which copies data from ``fin`` to ``fout`` in chunks of 4M.
    '''
    pid = os.fork()
    if pid == 0:
        os.close(pipe[0])
        checksum = hashlib.sha1()
        tail = ''
        length = 0
        while True:
            block = fin.read(COPY_BLOCK_SIZE)
            if len(block) == 0:
                break
                
            # update checksum
            checksum.update(block)
            
            # update tail
            if len(block) >= TAIL_LENGTH:
                tail = block[-TAIL_LENGTH:]
            else:
                keep_length = TAIL_LENGTH - len(block)
                tail = tail[0:keep_length] + block
                
            # update length
            length += len(block)
            
            # write block to output file
            if fout is not None:
                fout.write(block)
                
            # write block to pipe
            if use_pipe:
                bytes_written = os.write(pipe[1], block)
                
        # we're finished, close everything
        if fout != None:
            fout.close()
        if use_pipe:
            os.close(pipe[1])
            
        with open(report_path, 'w') as freport:
            report = dict()
            report['sha1'] = checksum.hexdigest()
            report['tail'] = tail
            report['length'] = length
            freport.write(yaml.dump(report))
            
        os._exit(0)
    else:
        name_for_pid[pid] = '[copy process]'
        add_proc_info(pid, name = "copy process to capture %s of PID %d" % (which, other_pid))
        proc_details[pid]['attach_to'] = [other_pid, which]
        log("Launched a copy process with PID " + str(pid) + " to capture " + which + " of PID " + str(other_pid) + ".")
        return pid


def launch(args, stdout_path = None, stderr_path = None, use_stdin = subprocess.PIPE):
    '''
    Launch a process.
    '''

    proc = subprocess.Popen(args,
        stdin = use_stdin,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        bufsize = COPY_BLOCK_SIZE,
        preexec_fn = os.setsid
    )
    name_for_pid[proc.pid] = os.path.basename(args[0])
    add_proc_info(proc.pid, args = args)
    log("Launched " + ' '.join(args) + " as PID " + str(proc.pid) + '.')

    pipe = os.pipe()
    
    stdout_sink = None
    if stdout_path != None:
        stdout_sink = open(stdout_path, 'w')
    report_path = temppath()
    pid = launch_copy_process(proc.stdout, stdout_sink, report_path, proc.pid, 'stdout', pipe, True)
    copy_process_reports[pid] = {
        'report_path': report_path,
        'stream': 'stdout',
        'pid': proc.pid,
        'args': args,
        'sink': os.path.basename(stdout_path) if stdout_path else None
    }

    os.close(pipe[1])
    
    stderr_sink = None
    if stderr_path != None:
        stderr_sink = open(stderr_path, 'w')
    report_path = temppath()
    pid = launch_copy_process(proc.stderr, stderr_sink, report_path, proc.pid, 'stderr', pipe, False)
    copy_process_reports[pid] = {
        'report_path': report_path,
        'stream': 'stderr',
        'pid': proc.pid,
        'args': args,
        'sink': os.path.basename(stderr_path) if stderr_path else None
    }
    
    return (proc, pipe[0])

def wait():
    '''
    Wait until all child processes have finished.
    '''
    
    # Seal all pipelines, i. e. add a consumer to the end
    global seal_these_pipeline_instances
    for p in seal_these_pipeline_instances:
        p.seal()
    seal_these_pipeline_instances = list()
        
    something_went_wrong = False
    while True:
        try:
            pid, exitcode = os.wait()
            log("%s (PID %d) has exited with exit code %d." % ((proc_details[pid]['name'] if 'name' in proc_details[pid] else proc_details[pid]['args']), pid, exitcode))
            proc_details[pid]['exit_code'] = exitcode
            try:
                if pid in copy_process_reports:
                    report_path = copy_process_reports[pid]['report_path']
                    if os.path.exists(report_path):
                        report = yaml.load(open(report_path, 'r'))
                        if copy_process_reports[pid]['sink']:
                            sha1_checksum_for_file_basename[copy_process_reports[pid]['sink']] = report['sha1']
                            proc_details[pid]['sink'] = copy_process_reports[pid]['sink']
                        copy_process_reports[pid].update(report)
                        del copy_process_reports[pid]['report_path']
                        os.unlink(report_path)
                        proc_details[pid].update(report)
            except Exception:
                # swallow possible exceptions because getting the checksum and
                # the tail is nice to have, but no cause for the entire process
                # to fail
                pass
        except TimeoutException:
            log("Timeout, killing all child processes now.")
            kill_all_child_processes()
        except OSError:
            # no more children running, we are done
            signal.alarm(0)
            log("Cancelling timeout (if there was one), all children have exited.")
            break
        else:
            if exitcode != 0:
                # Oops, something went wrong. See what happens and terminate
                # all child processes in a few seconds.
                something_went_wrong = True
                signal.signal(signal.SIGALRM, timeout_handler) 
                log("Terminating all children in %d seconds..." % SIGTERM_TIMEOUT)
                signal.alarm(SIGTERM_TIMEOUT)
                        
    if something_went_wrong:
        log("Pipeline crashed.")
        raise StandardError("Pipeline crashed.")
                    
class UnixPipeline(object):
    '''
    This class can be used to chain multiple processes together.
    '''
    def __init__(self):
        self.pipeline_procs = []
        self.use_stdin = subprocess.PIPE
        seal_these_pipeline_instances.append(self)

    def append(self, args, stdout_path = None, stderr_path = None):
        '''
        Append a process to the pipeline. If there already is a process in
        the pipeline, its ``stdout`` will be connected to the ``stdin`` of 
        the new process. File objects may be passed to capture ``stdout`` 
        and ``stderr``.
        '''
        proc, stdout_copy = launch(args, stdout_path, stderr_path, use_stdin = self.use_stdin)
        self.use_stdin = stdout_copy

        self.pipeline_procs.append(proc)

    def seal(self):
        pid = os.fork()
        if pid == 0:
            while True:
                block = os.read(self.use_stdin, COPY_BLOCK_SIZE)
                if len(block) == 0:
                    break
            os._exit(0)
        else:
            name_for_pid[pid] = 'null consumer'
            add_proc_info(pid, name = "null consumer for PID %d" % pid)
            log("Launched a null consumer with PID %d." % pid)
            