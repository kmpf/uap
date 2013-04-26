import sys
from abstract_step import *
import copy
import glob
import yaml

class Source(AbstractStep):

    connections = []
    connections.append('out/reads')

    def __init__(self, pipeline):
        super(Source, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        output_run_info = {}
        '''
        for key, sample_info in self.pipeline.all_samples.items():
            input_files = sample_info['files']
            output_run_info[key] = { 'output_files': { 'reads' : {} } }
            for path in input_files:
                output_run_info[key]['output_files']['reads'][path] = []
            if 'info' in sample_info:
                output_run_info[key]['info'] = copy.deepcopy(sample_info['info'])
        '''
        output_run_info['RIB0000740'] = dict()
        output_run_info['RIB0000740']['output_files'] = dict()
        output_run_info['RIB0000740']['output_files']['reads'] = dict()
        output_run_info['RIB0000740']['output_files']['reads']['/home/michael/Desktop/rnaseq-pipeline/in/Unaligned/Project_A/Sample_RIB0000740/RIB0000740_GATCAG_L001_R1_001-head.fastq.gz'] = []
        return output_run_info
