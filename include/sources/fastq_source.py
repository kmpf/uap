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

        self.add_option('indices', str, dict, optional = True)
        self.add_option('pattern', str)
        self.add_option('group', str)
        self.add_option('paired_end', bool)
        self.add_option('sample_id_prefix', str, optional = True)
                
    def setup_runs(self, input_run_info, connection_info):
        regex = re.compile(self.option('group'))

        output_run_info = {}

        for path in glob.glob(self.option('pattern')):
            match = regex.match(os.path.basename(path))
            sample_id_parts = []
            if self.option_set_in_config('sample_id_prefix'):
                sample_id_parts.append(self.option('sample_id_prefix'))
            sample_id_parts += list(match.groups())
            sample_id = '_'.join(sample_id_parts)
            if not sample_id in output_run_info:
                output_run_info[sample_id] = { 'output_files': { 'reads': {} }, 'info': { 'paired_end': self.option('paired_end') } }
            output_run_info[sample_id]['output_files']['reads'][path] = []

        if self.option_set_in_config('indices'):
            if type(self.option('indices')) == str:
                # read indices from CSV file or dictionary
                indices_path = self.option('indices')
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
                for sample_id, index in self.option('indices').items():
                    if sample_id in output_run_info:
                        if not 'index' in output_run_info[sample_id]['info']:
                            output_run_info[sample_id]['info']['index'] = index
                        else:
                            if index != output_run_info[sample_id]['info']['index']:
                                raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)

                        
        if self.option('paired_end'):
            # determine R1/R2 info for each input file: read_number
            for sample_name in output_run_info.keys():
                output_run_info[sample_name]['info']['read_number'] = {}
                for path in output_run_info[sample_name]['output_files']['reads'].keys():
                    isR1 = '_R1' in path or '-R1' in path
                    isR2 = '_R2' in path or '-R2' in path
                    if isR1 and isR2:
                        raise StandardError("Unable to determine read_numer, seems to be both R1 and R2: " + path)
                    if (not isR1) and (not isR2):
                        raise StandardError("Unable to determine read_numer, seems to be neither R1 nor R2: " + path)
                    output_run_info[sample_name]['info']['read_number'][os.path.basename(path)] = 'R1' if isR1 else 'R2'
                    
        return output_run_info
