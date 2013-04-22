import sys
sys.path.append('./include/steps')
import os
import signal
import subprocess
import tempfile

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

up_log = []

def log(message):
    up_log.append(message)
    
def get_log():
    return "\n".join(up_log)
    
def mkfifo(id):
    _, path = tempfile.mkstemp(id)
    os.close(_)
    os.unlink(path)
    os.mkfifo(path)
    return path

def kill_all_child_processes():
    for pid, name in name_for_pid.items():
        try:
            os.kill(pid, signal.SIGTERM)
            sys.stderr.write("Killed " + name + " (PID " + str(pid) + ").\n")
            sys.stderr.flush()
        except OSError:
            pass


def launch_copy_thread(fin, fout):
    pid = os.fork()
    if pid == 0:
        while True:
            block = fin.read(4096 * 1024)
            if len(block) == 0:
                break
            fout.write(block)
        fout.close()
        os._exit(0)
    else:
        name_for_pid[pid] = '[copy thread]'
        log("[up] Launched a copy thread with PID " + str(pid) + ".")


def launch(args, stdout = None, stderr = None, use_stdin = subprocess.PIPE):
    proc = subprocess.Popen(args,
        stdin = use_stdin,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        bufsize = 4096 * 1024,
        preexec_fn = os.setsid
    )
    name_for_pid[proc.pid] = os.path.basename(args[0])
    log("[up] Launched " + ' '.join(args) + " as PID " + str(proc.pid) + '.')

    if stdout != None:
        launch_copy_thread(proc.stdout, stdout)
    if stderr != None:
        launch_copy_thread(proc.stderr, stderr)
    return proc

def wait():
    # wait until all child processes have finished
    while True:
        try:
            pid, exitcode = os.wait()
            log("[up] Child " + str(pid) + " has exited with exit code " + str(exitcode) + ".")
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
                message += "Full pipeline log:\n\n" + get_log() + "\n"
                raise StandardError(message)
        else:
            if pid in upstream_procs:
                for upstream_proc in upstream_procs[pid]:
                    ok_to_fail.append(upstream_proc.pid)
                    try:
                        upstream_proc.terminate()
                    except OSError:
                        pass
                    
def create_pipeline():
    return UnixPipeline()
                    
class UnixPipeline(object):
    def __init__(self):
        self.pipeline_procs = []
        self.use_stdin = None

    def append(self, args, stdout = None, stderr = None):
        if len(self.pipeline_procs) > 0:
            self.use_stdin = self.pipeline_procs[-1].stdout
            
        proc = launch(args, stdout, stderr, use_stdin = self.use_stdin)

        self.pipeline_procs.append(proc)
        upstream_procs[proc.pid] = self.pipeline_procs[0:-1]
