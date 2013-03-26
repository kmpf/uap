import sys
from abstract_step import *
import pipeline
import subprocess
import yaml

class Cutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)

    def setup_runs(self, complete_input_run_info):
        output_run_info = {}
        for input_run_id, input_run_info in complete_input_run_info.items():
            for in_path in sorted(input_run_info['output_files']['reads'].keys()):
                # determine which read this is (R1 or R2)
                which = None
                if '_R1_' in in_path:
                    which = 'R1'
                elif '_R2_' in in_path:
                    which = 'R2'
                else:
                    raise StandardError("Expected input files with _R1_ or _R2_.")

                output_run_id = input_run_id + '-' + which

                if not output_run_id in output_run_info:
                    output_run_info[output_run_id] = {
                        'output_files': {},
                        'info': {
                            'read': which
                        }
                    }

                for t in [('reads', input_run_id + '-cutadapt-' + which + '.fastq.gz'),
                        ('log', input_run_id + '-cutadapt-' + which + '-log.txt')]:
                    pathkey = t[0]
                    path = t[1]
                    if not pathkey in output_run_info[output_run_id]['output_files']:
                        output_run_info[output_run_id]['output_files'][pathkey] = {}
                    if not path in output_run_info[output_run_id]['output_files'][pathkey]:
                        output_run_info[output_run_id]['output_files'][pathkey][path] = []
                    output_run_info[output_run_id]['output_files'][pathkey][path].append(in_path)

        return output_run_info

    def execute(self, run_id, run_info):
        if len(run_info['output_files']['reads']) != 1:
            raise StandardError("Expected a single output file.")

        basename = os.path.basename(run_info.keys()[0])
        adapter = ''
        if run_info['info']['read'] == 'R1':
            adapter = self.options['adapter-R1']
        elif run_info['info']['read'] == 'R2':
            adapter = self.options['adapter-R2']

        # TODO: replace (__INDEX__) in adapter
        if '(__INDEX__)' in adapter:
            # TODO: this is weird, we need something more general
            sample_info = self.pipeline.all_samples[run_id[0:-3]]
            index = sample_info['lanes'].values()[0]['Index']
            adapter = adapter.replace('(__INDEX__)', index)

        pigz1 = [self.pipeline.config['tools']['pigz']['path'], '-d', '-c']
        pigz1.extend(*sorted(run_info['output_files']['reads'].values()))

        cutadapt = [self.pipeline.config['tools']['cutadapt']['path'], '-a',
            adapter, '-']

        pigz2 = [self.pipeline.config['tools']['pigz']['path'],
            '--blocksize', '4096', '--processes', '3', '-c']

        procs = []
        procs_pid = []
        use_stdin = None
        for args in [pigz1, cutadapt, pigz2]:
            if len(procs) > 0:
                use_stdin = procs[-1].stdout
            proc = subprocess.Popen(args,
                    stdout = subprocess.PIPE,
                    bufsize = 4096 * 1024,
                    stdin = use_stdin,
                    stderr = subprocess.PIPE,
                    preexec_fn = os.setsid
                    )
            procs.append(proc)
            procs_pid.append(proc.pid)

        copy_streams = []
        copy_streams.append((procs[-1].stdout, open(run_info['output_files']['reads'].keys()[0], 'w')))
        copy_streams.append((procs[-2].stderr, open(run_info['output_files']['log'].keys()[0], 'w')))

        # copy data in threads
        for info in copy_streams:
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
                procs_pid.append(pid)

        while procs_pid:
            try:
                pid, exitcode = os.wait()
            except OSError:
                break
            if exitcode != 0:
                print(pid)
                print(exitcode)
                raise StandardError("PIPELINE CRASHED, OH MY.")
            procs_pid.remove(pid)
