import sys
sys.path.append('./include/steps')
import os
import subprocess

class UnixPipeline(object):
    def __init__(self):
        self.procs = []

        # the stdout of the last process
        self.use_stdin = None
        self.copy_streams = []
        self.upstream_procs = {}
        self.ok_to_fail = []

    def append(self, args, stdout = None, stderr = None):
        if len(self.procs) > 0:
            self.use_stdin = self.procs[-1].stdout
        sys.stderr.write("[up] Launching " + ' '.join(args) + " ... ")
        proc = subprocess.Popen(args,
            stdout = subprocess.PIPE,
            bufsize = 4096 * 1024,
            stdin = self.use_stdin,
            stderr = subprocess.PIPE,
            preexec_fn = os.setsid
        )
        sys.stderr.write("launched as PID " + str(proc.pid) + "\n")
        self.procs.append(proc)
        self.upstream_procs[proc.pid] = self.procs[0:-1]

        if stdout != None:
            self.copy_streams.append((proc.stdout, stdout))
        if stderr != None:
            self.copy_streams.append((proc.stderr, stderr))

    def run(self):
        # set up threads for writing data
        for info in self.copy_streams:
            pid = os.fork()
            if pid == 0:
                while True:
                    block = info[0].read(4096 * 1024)
                    if len(block) == 0:
                        break
                    info[1].write(block)
                info[1].close()
                os._exit(0)
            else:
                sys.stderr.write("[up] Launched a copy thread with PID " + str(pid) + ".\n")

        # wait until all child processes have finished
        while True:
            try:
                pid, exitcode = os.wait()
                sys.stderr.write("[up] Child " + str(pid) + " has exited with exit code " + str(exitcode) + ".\n")
            except OSError:
                break
            if exitcode != 0:
                if not pid in self.ok_to_fail:
                    sys.stdout.flush()
                    sys.stderr.flush()
                    raise StandardError("PIPELINE CRASHED, OH MY.")
            else:
                if pid in self.upstream_procs:
                    for upstream_proc in self.upstream_procs[pid]:
                        self.ok_to_fail.append(upstream_proc.pid)
                        try:
                            upstream_proc.terminate()
                        except OSError:
                            pass
