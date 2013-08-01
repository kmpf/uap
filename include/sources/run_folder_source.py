import sys
from abstract_step import *
import csv
import glob
import os
import pipeline
import yaml

class RunFolderSource(AbstractSourceStep):
    def __init__(self, pipeline):
        super(RunolderSource, self).__init__(pipeline)
        
        self.add_connection('out/reads')
        
        self.add_option('path', str)
        self.add_option('paired_end', bool)
        
    def setup_runs(self, input_run_info, connection_info):
        path = self.get_option('path')

        output_run_info = {}

        if not os.path.exists(path):
            raise ConfigurationException("Source path does not exist: " + path)

        # find all samples
        for sample_path in glob.glob(os.path.join(path, 'Unaligned', 'Project_*', 'Sample_*')):
            sample_name = os.path.basename(sample_path).replace('Sample_', '')
            if sample_name in output_run_info:
                raise ConfigurationException("Duplicate sample: " + sample_name)
            
            if not sample_name in output_run_info:
                output_run_info[sample_name] = { 'output_files': { 'reads': {} }, 'info': { 'paired_end': self.option('paired_end') } }
            
            for path in sorted(glob.glob(os.path.join(sample_path, '*.fastq.gz'))):
                output_run_info[sample_name]['output_files']['reads'][path] = []

            # read sample sheets
            sample_sheet_path = os.path.join(sample_path, 'SampleSheet.csv')
            reader = csv.DictReader(open(sample_sheet_path))
            for row in reader:
                sample_id = row['SampleID']
                index = row['Index']
                if not 'index' in output_run_info[sample_id]['info']:
                    output_run_info[sample_id]['info']['index'] = index
                else:
                    if index != output_run_info[sample_id]['info']['index']:
                        raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)

            if self.option('paired_end') == True:
                # determine R1/R2 info for each input file: read_number
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
