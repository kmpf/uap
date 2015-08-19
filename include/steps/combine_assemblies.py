import sys
import re
import yaml

from ..abstract_step import *
from .. import pipeline
from .. import process_pool

class Combine_Assemblies(AbstractStep):
    
    def __init__(self, pipeline):
        super(Combine_Assemblies, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/features_tophat', constraints = {'total_files': 1} )
        self.add_connection('in/features_segemehl', constraints = {'total_files': 1} )
        
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        self.add_connection('out/features')
        
        self.require_tool('Rscript')

        self.add_option('run_id', str)
        self.add_option('mode', str, choices = ['no_transcript_overlap'])
        self.add_option('path_rscript', str)

    def declare_runs(self):
        features_segemehl_path = self.get_single_input_file_for_connection('in/features_segemehl')
        features_tophat_path = self.get_single_input_file_for_connection('in/features_tophat')
        run_id = self.get_option('run_id')
        with self.declare_run(run_id) as run:
            run.add_private_info('features_segemehl_path', features_segemehl_path)
            run.add_private_info('features_tophat_path', features_tophat_path)
            run.add_output_file('features', '%s.gtf' % run_id, [features_segemehl_path, features_tophat_path])
            run.add_output_file('log_stderr', '%s-log_stderr.txt' % run_id, [features_segemehl_path, features_tophat_path])
            run.add_output_file('log_stdout', '%s-log_stdout.txt' % run_id, [features_segemehl_path, features_tophat_path])
            run.new_exec_group()
            
    def execute(self, run_id, run):
        
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                
                ## Append the input files, which are all input file linked to the output file
                item_A_path = run.get_private_info('features_segemehl_path')
                item_B_path = run.get_private_info('features_tophat_path')

                ## Obtain the output files
                gtf_out_path = run.get_single_output_file_for_annotation('features')

                ## Call Rscript
                rscript = [self.get_tool('Rscript'), self.get_option('path_rscript')]

                ## Append input files
                rscript.extend (['--inputA', item_A_path])
                rscript.extend (['--inputB', item_B_path])

                ## Append the output files
                rscript.extend (['--output', gtf_out_path])

                ## Append prefixes
                rscript.extend (['--prefixA', 'segemehl'])
                rscript.extend (['--prefixB', 'tophat'])

                ## Append mode of combining the GTF files
                rscript.extend (['--mode', 'no_transcript_overlap'])

                log_stderr = run.get_single_output_file_for_annotation('log_stderr')
                log_stdout = run.get_single_output_file_for_annotation('log_stdout')
                    
                pipeline.append (rscript, stderr_path = log_stderr, stdout_path = log_stdout)


