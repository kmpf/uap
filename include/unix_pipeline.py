import sys
sys.path.append('./include/steps')
import os
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


def mkfifo(id):
    _, path = tempfile.mkstemp(id)
    os.close(_)
    os.unlink(path)
    os.mkfifo(path)
    return path


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
        sys.stderr.write("[up] Launched a copy thread with PID " + str(pid) + ".\n")


def launch(args, stdout = None, stderr = None):
    sys.stderr.write("[up] Launching " + ' '.join(args) + " ... ")
    proc = subprocess.Popen(args,
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        bufsize = 4096 * 1024,
        preexec_fn = os.setsid
    )
    sys.stderr.write("launched as PID " + str(proc.pid) + "\n")

    if stdout != None:
        launch_copy_thread(proc.stdout, stdout)
    if stderr != None:
        launch_copy_thread(proc.stderr, stderr)

def wait():
    # wait until all child processes have finished
    while True:
        try:
            pid, exitcode = os.wait()
            sys.stderr.write("[up] Child " + str(pid) + " has exited with exit code " + str(exitcode) + ".\n")
        except OSError:
            # no more children running, we are done
            break
        if exitcode != 0:
            if not pid in ok_to_fail:
                # oops, something went wrong
                sys.stdout.flush()
                sys.stderr.flush()
                raise StandardError("PIPELINE CRASHED, OH MY.")
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
        sys.stderr.write("[up] Launching " + ' '.join(args) + " ... ")
        proc = subprocess.Popen(args,
            stdout = subprocess.PIPE,
            bufsize = 4096 * 1024,
            stdin = self.use_stdin,
            stderr = subprocess.PIPE,
            preexec_fn = os.setsid
        )
        sys.stderr.write("launched as PID " + str(proc.pid) + "\n")
        self.pipeline_procs.append(proc)
        upstream_procs[proc.pid] = self.pipeline_procs[0:-1]
        
        if stdout != None:
            launch_copy_thread(proc.stdout, stdout)
        if stderr != None:
            launch_copy_thread(proc.stderr, stderr)
           