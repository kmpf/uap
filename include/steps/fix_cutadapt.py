import sys
from abstract_step import *

class FixCutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)

    def setup_runs(self, input_run_info_complete):
        output_run_info = {}
        for input_run_id, input_run_info in input_run_info_complete.items():
            new_key = input_run_id.replace('-R1', '').replace('-R2', '')
            if not new_key in output_run_info:
                output_run_info[new_key] = { 'output_files': { 'reads': {} } }
            for in_path in sorted(input_run_info['output_files']['reads'].keys()):
                for _ in ['R1', 'R2']:
                    k2 = new_key + '-fixed-' + _ + '.fastq.gz'
                    if not k2 in output_run_info[new_key]['output_files']['reads']:
                        output_run_info[new_key]['output_files']['reads'][k2] = []
                    output_run_info[new_key]['output_files']['reads'][k2].append(in_path)
        return output_run_info