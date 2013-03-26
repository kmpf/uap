import sys
from abstract_step import *
import pipeline
import pooryorick_pipeline

class Head(AbstractStep):
    def __init__(self, pipeline):
        super(Head, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        output_run_info = {}
        for key, input_files in input_run_info.items():
            for fn in input_files:
                k = key
                if not k in output_run_info:
                    output_run_info[k] = {}
                destination_file_name = fn.replace('.fastq.gz', '') + '-head.fastq.gz'
                if not destination_file_name in output_run_info[k]:
                    output_run_info[k][destination_file_name] = []
                output_run_info[k][destination_file_name].append(fn)
        return output_run_info

    def execute(self, run_id, run_info):
        for outpath, inpaths in run_info.items():
            if len(inpaths) != 1:
                raise StandardError("Expected one input file per output file.")
            inpath = inpaths[0]

            pigz1 = [self.pipeline.config['tools']['pigz']['path'], '-d', '-c', inpath]

            head = ['head', '-n', '1000']

            pigz2 = [self.pipeline.config['tools']['pigz']['path'],
                '--blocksize', '4096', '--processes', '3', '-c']

            with open(outpath, 'w') as fout:
                fd, pids = pooryorick_pipeline.pipeline(pigz1, head, pigz2)
                while True:
                    block = os.read(fd, 4096 * 1024)
                    if len(block) == 0:
                        break
                    fout.write(block)