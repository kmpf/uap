'''
This module can be used to launch child processes and wait for them.
Processes may either run on their own or pipelines can be built with them.
'''

import sys
sys.path.append('./include/steps')
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

# a list of PIDs which are ok to fail because a child has terminated
# for example, if we do:
# $ cat [file] | head -n 10
# ...head would stop after ten lines. cat would run on forever if we
# wouldn't terminate it, so we do. But before that, we must remember
# that it's OK for cat to fail. Right?
# Values are integers (PIDs).
ok_to_fail = []

# a dict of upstream procs for any proc in a pipeline, values are Popen objects
# A => B => C
# upstream_procs:
#   B: [A]
#   C: [B, A]
upstream_procs = {}

# for better error messages, this dict holds the process names for PIDs
name_for_pid = {}

sha1_checksum_for_file = {}

up_log = []
copy_thread_reports = {}
pipeline_instances = []

def clear():
    ok_to_fail = list()
    upstream_proc = dict()
    name_for_pid = dict()
    sha1_checksum_for_file = dict()
    up_log = list()
    copy_thread_reports = dict()
    pipeline_instances = list()

def log(message):
    '''
    Append a message to the pipeline log.
    '''
    up_log.append("[up] [%s] %s" % (datetime.datetime.now(), message))
    
def get_log():
    '''
    Return the log as a dictionary.
    '''
    log = dict()
    log['messages'] = up_log
    log['streams'] = copy_thread_reports.values()
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
            freport.write(yaml.dump(report))
            
        os._exit(0)
    else:
        name_for_pid[pid] = '[copy process]'
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
    log("Launched " + ' '.join(args) + " as PID " + str(proc.pid) + '.')

    pipe = os.pipe()
    
    stdout_sink = None
    if stdout_path != None:
        stdout_sink = open(stdout_path, 'w')
    report_path = temppath()
    pid = launch_copy_process(proc.stdout, stdout_sink, report_path, proc.pid, 'stdout', pipe, True)
    copy_thread_reports[pid] = {
        'report_path': report_path,
        'stream': 'stdout',
        'pid': proc.pid,
        'args': args,
        'sink': stdout_path
    }

    os.close(pipe[1])
    
    stderr_sink = None
    if stderr_path != None:
        stderr_sink = open(stderr_path, 'w')
    report_path = temppath()
    pid = launch_copy_process(proc.stderr, stderr_sink, report_path, proc.pid, 'stderr', pipe, False)
    copy_thread_reports[pid] = {
        'report_path': report_path,
        'stream': 'stderr',
        'pid': proc.pid,
        'args': args,
        'sink': stderr_path
    }
    
    return (proc, pipe[0])

def wait():
    '''
    Wait until all child processes have finished.
    '''
    
    # Seal all pipelines, i. e. add a consumer to the end
    for p in pipeline_instances:
        p.seal()
        
    something_went_wrong = None
    while True:
        try:
            pid, exitcode = os.wait()
            log("Child " + str(pid) + " has exited with exit code " + str(exitcode) + ".")
            try:
                if pid in copy_thread_reports:
                    report_path = copy_thread_reports[pid]['report_path']
                    if os.path.exists(report_path):
                        report = yaml.load(open(report_path, 'r'))
                        if copy_thread_reports[pid]['sink']:
                            sha1_checksum_for_file[copy_thread_reports[pid]['sink']] = report['sha1']
                        copy_thread_reports[pid].update(report)
                        del copy_thread_reports[pid]['report_path']
                        os.unlink(report_path)
            except:
                # swallow possible exceptions because getting the checksum and
                # the tail is nice to have, but no cause for the entire process
                # to fail
                pass
        except OSError:
            # no more children running, we are done
            break
        if exitcode != 0:
            if not pid in ok_to_fail:
                # oops, something went wrong
                sys.stdout.flush()
                sys.stderr.flush()
                message = "Pipeline crashed, oh my.\n"
                job_name = 'Process with PID ' + str(pid)
                if pid in name_for_pid:
                    job_name = name_for_pid[pid] + ' (PID ' + str(pid) + ')'
                message += job_name + ' has crashed with exit code ' + str(exitcode) + '.\n'
                message += "Full pipeline log:\n\n%s\n" % yaml.dump(get_log(), default_flow_style = False)
                something_went_wrong = message
            else:
                if pid in upstream_procs:
                    for upstream_proc in upstream_procs[pid]:
                        ok_to_fail.append(upstream_proc.pid)
                        try:
                            upstream_proc.terminate()
                        except OSError:
                            pass
                        
    if something_went_wrong != None:
        sys.stderr.write(something_went_wrong)
        raise StandardError("Pipeline crashed.")
                    
class UnixPipeline(object):
    '''
    This class can be used to chain multiple processes together.
    '''
    def __init__(self):
        self.pipeline_procs = []
        self.use_stdin = subprocess.PIPE
        pipeline_instances.append(self)

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
        upstream_procs[proc.pid] = self.pipeline_procs[0:-1]

    def seal(self):
        if os.fork() == 0:
            while True:
                block = os.read(self.use_stdin, COPY_BLOCK_SIZE)
                if len(block) == 0:
                    break
            os._exit(0)