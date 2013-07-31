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

    def declare_runs(self):
        for run_id, input_paths in self.run_ids_and_input_files_for_connection('in/*'):
            with self.declare_run(run_id) as run:
                for in_path in input_paths:
                    annotation = self.annotation_for_input_file(in_path)
                    self.add_connection('out/%s' % annotation)
                    run.add_output_file(annotation, misc.append_suffix_to_path(os.path.basename(in_path), 'head'), [in_path])
    
    def execute(self, run_id, run):
        # process one file at a time
        for tag, output_file_info in run.output_files().items():
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
