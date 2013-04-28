import sys
from abstract_step import *
import copy
import csv
import glob
import os
import re
import yaml

class FastqSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(FastqSource, self).__init__(pipeline)
        
        self.add_connection('out/reads')

    def setup_runs(self, input_run_info):
        regex = re.compile(self.options['group'])

        output_run_info = {}

        if not 'paired_end' in self.options:
            raise StandardError("missing paired_end key in source")

        for path in glob.glob(self.options['pattern']):
            match = regex.match(os.path.basename(path))
            sample_id_parts = []
            if 'sample_id_prefix' in self.options:
                sample_id_parts.append(self.options['sample_id_prefix'])
            sample_id_parts += list(match.groups())
            sample_id = '_'.join(sample_id_parts)
            if not sample_id in output_run_info:
                output_run_info[sample_id] = { 'output_files': { 'reads': {} }, 'info': { 'paired_end': self.options['paired_end'] } }
            output_run_info[sample_id]['output_files']['reads'][path] = []

        if self.options['indices'].__class__ == str:
            # read indices from CSV file or dictionary
            indices_path = self.options['indices']
            reader = csv.DictReader(open(indices_path))
            for row in reader:
                sample_id = row['SampleID']
                if sample_id in output_run_info:
                    index = row['Index']
                    if not 'index' in output_run_info[sample_id]['info']:
                        output_run_info[sample_id]['info']['index'] = index
                    else:
                        if index != output_run_info[sample_id]['info']['index']:
                            raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)
        else:
            # indices are defined in the configuration
            for sample_id, index in self.options['indices'].items():
                if sample_id in output_run_info:
                    if not 'index' in output_run_info[sample_id]['info']:
                        output_run_info[sample_id]['info']['index'] = index
                    else:
                        if index != output_run_info[sample_id]['info']['index']:
                            raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)

        if self.options['paired_end'] == True:
            # determine R1/R2 info for each input file: read_number
            for sample_name in output_run_info.keys():
                output_run_info[sample_name]['info']['read_number'] = {}
                for path in output_run_info[sample_name]['output_files']['reads'].keys():
                    isR1 = '_R1' in path
                    isR2 = '_R2' in path
                    if isR1 and isR2:
                        raise StandardError("Unable to determine read_numer, seems to be both R1 and R2: " + path)
                    if (not isR1) and (not isR2):
                        raise StandardError("Unable to determine read_numer, seems to be neither R1 nor R2: " + path)
                    output_run_info[sample_name]['info']['read_number'][os.path.basename(path)] = 'R1' if isR1 else 'R2'
                    
        return output_run_info
