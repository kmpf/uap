import sys
from abstract_step import *

class FixCutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        output_run_info = {}
        for key, input_files in input_run_info.items():
            new_key = key.replace('-R1', '').replace('-R2', '')
            if not new_key in output_run_info:
                output_run_info[new_key] = {}
            for fn in input_files:
                for _ in ['R1', 'R2']:
                    k2 = new_key + '-fixed-' + _ + '.fastq.gz'
                    if not k2 in output_run_info[new_key]:
                        output_run_info[new_key][k2] = []
                    output_run_info[new_key][k2].append(fn)
        return output_run_info