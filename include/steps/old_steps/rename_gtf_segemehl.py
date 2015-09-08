import sys
from abstract_step import *
import process_pool
import re
import yaml


class Rename_GTF_Segemehl(AbstractStep):

    def __init__(self, pipeline):
        super(Rename_GTF_Segemehl, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/features')
        self.add_connection('out/features_segemehl')
        
        self.require_tool('cat')

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/features'):
            with self.declare_run(run_id) as run:
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file., but got this %s" % input_paths)
                run.add_private_info('in-features', input_paths[0])
                run.add_output_file('features_segemehl', '%s-segemehl.gtf' % run_id, input_paths)
                run.new_exec_group()
                
    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                features_path = run.get_private_info('in-features')
                cat = [self.get_tool('cat'), features_path]
                               
                pipeline.append(cat, stdout_path = run.get_single_output_file_for_annotation('features_segemehl'))
