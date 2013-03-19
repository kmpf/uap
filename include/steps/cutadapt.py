import sys
from abstract_step import *
import yaml

class Cutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)

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
                elif '_R2_' in fn:
                    k = key + '-R2'
                    if not k in output_run_info:
                        output_run_info[k] = {}
                    destination_file_name = key + '-cutadapt-R2.fastq.gz'
                    if not destination_file_name in output_run_info[k]:
                        output_run_info[k][destination_file_name] = []
                    output_run_info[k][destination_file_name].append(fn)
                else:
                    raise StandardError("Expected input files with _R1_ or _R2_.")
        return output_run_info

    def command_line(self, key):
        command = ''
        command = "\"" + self.options['path'] + "\""
