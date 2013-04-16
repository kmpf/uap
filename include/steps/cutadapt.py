import sys
from abstract_step import *
import pipeline
import re
import subprocess
import unix_pipeline
import yaml

class Cutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)
        self.set_cores(6)

    def setup_runs(self, complete_input_run_info):
        # make sure tools are available
        self.tool('pigz')
        self.tool('cutadapt')

        output_run_info = {}
        for input_run_id, input_run_info in complete_input_run_info.items():
            for in_path in sorted(input_run_info['output_files']['reads'].keys()):
                suffix = ''
                which = None
                if input_run_info['info']['paired_end'] == True:
                    which = input_run_info['info']['read_number'][os.path.basename(in_path)]
                    if not which in ['R1', 'R2']:
                        raise StandardError("Expected R1 and R2 input files, but got this: " + in_path)
                    suffix = '-' + which

                output_run_id = input_run_id + suffix

                if not output_run_id in output_run_info:
                    output_run_info[output_run_id] = {
                        'output_files': {},
                        'info': {}
                    }
                    if input_run_info['info']['paired_end'] == True:
                        output_run_info[output_run_id]['info']['read_number'] = which

                # find adapter
                adapter = self.options['adapter' + suffix]

                # insert correct index if necessary
                if '((INDEX))' in adapter:
                    index = input_run_info['info']['index']
                    adapter = adapter.replace('((INDEX))', index)

                # make sure the adapter is looking good
                if re.search('^[ACGT]+$', adapter) == None:
                    raise StandardError("Unable to come up with a legit-looking adapter: " + adapter)
                output_run_info[output_run_id]['info']['adapter'] = adapter

                for t in [('reads', input_run_id + '-cutadapt' + suffix + '.fastq.gz'),
                        ('log', input_run_id + '-cutadapt' + suffix + '-log.txt')]:
                    pathkey = t[0]
                    path = t[1]
                    if not pathkey in output_run_info[output_run_id]['output_files']:
                        output_run_info[output_run_id]['output_files'][pathkey] = {}
                    if not path in output_run_info[output_run_id]['output_files'][pathkey]:
                        output_run_info[output_run_id]['output_files'][pathkey][path] = []
                    output_run_info[output_run_id]['output_files'][pathkey][path].append(in_path)

        return output_run_info

    def execute(self, run_id, run_info):
        # basic sanity check
        if len(run_info['output_files']['reads']) != 1:
            raise StandardError("Expected a single output file.")

        # set up processes
        cat4m = [self.tool('cat4m')]
        cat4m.extend(*sorted(run_info['output_files']['reads'].values()))

        pigz1 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '1', '-d', '-c']

        cutadapt = [self.tool('cutadapt'), '-a', run_info['info']['adapter'], '-']

        pigz2 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '3', '-c']

        # create the pipeline and run it
        up = unix_pipeline.UnixPipeline()
        up.append(cat4m)
        up.append(pigz1)
        up.append(cutadapt, stderr = open(run_info['output_files']['log'].keys()[0], 'w'))
        up.append(pigz2, stdout = open(run_info['output_files']['reads'].keys()[0], 'w'))

        up.run()
