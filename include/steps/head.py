import sys
from abstract_step import *
import copy
import pipeline
import re
import process_pool
import yaml

class Head(AbstractStep):
    
    def __init__(self, pipeline):
        super(Head, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/*')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('head')
        
        self.add_option('lines', int, default = 1000)

    def setup_runs(self, complete_input_run_info, connection_info):
        
        output_run_info = {}
        
        for step_id, step_info in complete_input_run_info.items():
            for run_id, run_info in step_info.items():
                output_run_info[run_id] = dict()
                output_run_info[run_id]['output_files'] = dict()
                for tag, output_files in run_info['output_files'].items():
                    self.add_connection('out/' + tag)
                    output_run_info[run_id]['output_files'][tag] = dict()
                    for output_path, input_paths in output_files.items():
                        output_run_info[run_id]['output_files'][tag][misc.append_suffix_to_path(output_path, 'head')] = [output_path]
        
        return output_run_info
    
    
    def execute(self, run_id, run_info):
        # process one file at a time
        for tag, output_file_info in run_info['output_files'].items():
            for output_path, input_paths in output_file_info.items():
                with process_pool.ProcessPool(self) as pool:
                    with pool.Pipeline(pool) as pipeline:
                        cat4m = [self.tool('cat4m'), input_paths[0]]
                        pigz1 = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                        head = [self.tool('head'), '-n', str(self.option('lines'))]
                        pigz2 = [self.tool('pigz'), '--processes', '2', '--blocksize', '4096', '--stdout']
                
                        if output_path[-3:] == '.gz':
                            pipeline.append(cat4m)
                            pipeline.append(pigz1)
                            pipeline.append(head)
                            pipeline.append(pigz2, stdout_path = output_path)
                        else:
                            pipeline.append(cat4m)
                            pipeline.append(head, stdout_path = output_path)
