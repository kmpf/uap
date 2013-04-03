import sys
from abstract_step import *
import glob
import yaml

class Source(AbstractStep):
    def __init__(self, pipeline):
        super(Source, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        output_run_info = {}
        for key, sample_info in self.pipeline.all_samples.items():
            input_files = sample_info['files']
            output_run_info[key] = { 'output_files': { 'reads' : {} } }
            for path in input_files:
                output_run_info[key]['output_files']['reads'][path] = []
        return output_run_info

    def get_run_ids(self):
        # because this is a source, we don't need to actually perform any
        # steps - return an empty list
        return []

    def command_line(self, key):
        return None
