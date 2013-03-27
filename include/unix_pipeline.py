import sys
sys.path.append('./include/steps')
import os
import subprocess

class UnixPipeline(object):
    def __init__(self):
        self.procs = []
        self.procs_pid = []
        self.use_stdin = None
        self.copy_streams = []
        self.upstream_procs = {}

    def append(self, args, stdout = None, stderr = None):
        if len(self.procs) > 0:
            self.use_stdin = self.procs[-1].stdout
        proc = subprocess.Popen(args,
            stdout = subprocess.PIPE,
            bufsize = 4096 * 1024,
            stdin = self.use_stdin,
            stderr = subprocess.PIPE,
            preexec_fn = os.setsid
        )
        self.procs.append(proc)
        self.procs_pid.append(proc.pid)
        self.upstream_procs[proc.pid] = self.procs[0:-1]

        if stdout != None:
            self.copy_streams.append((proc.stdout, stdout))
        if stderr != None:
            self.copy_streams.append((proc.stderr, stderr))

    def run(self):
        # set up threads for writing data
        for info in self.copy_streams:
            pid = os.fork()
            if not pid:
                while True:
                    block = info[0].read(4096 * 1024)
                    if len(block) == 0:
                        break
                    info[1].write(block)
                info[1].close()
                os._exit(0)
            else:
                self.procs_pid.append(pid)

        # wait until all processes have finished
        while self.procs_pid:
            try:
                pid, exitcode = os.wait()
            except OSError:
                break
            if exitcode != 0:
                print(pid)
                print(exitcode)
                raise StandardError("PIPELINE CRASHED, OH MY.")
            else:
                if pid in self.upstream_procs:
                    for upstream_proc in self.upstream_procs[pid]:
                        upstream_proc.kill()
            self.procs_pid.remove(pid)
