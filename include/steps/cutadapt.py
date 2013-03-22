import sys
from abstract_step import *
import pipeline
import pooryorick_pipeline
import yaml

class Cutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)

        # for every output file, remember which read it was (R1 or R2)
        self.options_for_output_file = {}

    def setup_runs(self, input_run_info):
        output_run_info = {}
        for key, input_files in input_run_info.items():
            for fn in input_files:
                if '_R1_' in fn:
                    k = key + '-R1'
                    if not k in output_run_info:
                        output_run_info[k] = {}
                    destination_file_name = key + '-cutadapt-R1.fastq.gz'
                    if not destination_file_name in output_run_info[k]:
                        output_run_info[k][destination_file_name] = []
                    output_run_info[k][destination_file_name].append(fn)
                    self.options_for_output_file[destination_file_name] = {'read': 'R1'}
                elif '_R2_' in fn:
                    k = key + '-R2'
                    if not k in output_run_info:
                        output_run_info[k] = {}
                    destination_file_name = key + '-cutadapt-R2.fastq.gz'
                    if not destination_file_name in output_run_info[k]:
                        output_run_info[k][destination_file_name] = []
                    output_run_info[k][destination_file_name].append(fn)
                    self.options_for_output_file[destination_file_name] = {'read': 'R2'}
                else:
                    raise StandardError("Expected input files with _R1_ or _R2_.")
        return output_run_info

    def execute(self, run_id, run_info):
        if len(run_info) != 1:
            raise StandardError("Expected a single output file.")

        basename = os.path.basename(run_info.keys()[0])
        options = self.options_for_output_file[basename]
        adapter = ''
        if options['read'] == 'R1':
            adapter = self.options['adapter-R1']
        elif options['read'] == 'R2':
            adapter = self.options['adapter-R2']

        # TODO: replace ((INDEX)) in adapter

        pigz1 = [self.pipeline.config['tools']['pigz']['path'], '-d', '-c']
        pigz1.extend(run_info.values()[0])

        cutadapt = [self.pipeline.config['tools']['cutadapt']['path'], '-a',
            adapter, '-']

        pigz2 = [self.pipeline.config['tools']['pigz']['path'],
            '--blocksize', '4096', '--processes', '3', '-c']

        with open(run_info.keys()[0], 'w') as fout:
            fd, pids = pooryorick_pipeline.pipeline(pigz1, cutadapt, pigz2)
            while True:
                block = os.read(fd, 4096 * 1024)
                if len(block) == 0:
                    break
                fout.write(block)
        s = pooryorick_pipeline.close(pids)
        if max(s) > 0:
            raise pooryorick_pipeline.PipelineError("return status:  %s" % s)
