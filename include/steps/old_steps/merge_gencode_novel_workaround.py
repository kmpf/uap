import sys
from abstract_step import *
import process_pool
import yaml


class  Merge_Gencode_Novel_Workaround(AbstractStep):

    def __init__(self, pipeline):
        super(Merge_Gencode_Novel_Workaround, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/features')
        self.add_connection('out/features')
        self.add_connection('out/log_stderr')

#        self.add_option('run_id', str)
        self.add_option('reference', str )
#        self.add_option('remove_gencode', bool,   default=False )
#        self.add_option('remove_unstranded', bool, default=False)

 
        self.require_tool('cat')
        self.require_tool('sort')
        
    def declare_runs(self):
     
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/features'):
            with self.declare_run(run_id) as run:
                run.new_exec_group()
                run.add_private_info('in-features', input_paths)
                run.add_output_file('features', '%s-full.gtf' % run_id, input_paths)
                run.add_output_file('log_stderr', '%s-log_stderr.txt' % run_id, input_paths)

        
        


        



                    
    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                features_path = run.get_private_info('in-features')

#                features_path.append( str(self.get_option('reference')))
                ref =self.get_option('reference')
                cat = [self.get_tool('cat'), features_path[0], ref ]

#                sort  = [self.get_tool('sort'),  '-k 1,1 -k 4g,4 -k 5g,5']
                sort  = [self.get_tool('sort'),  '-k','1,1','-k', '4g,4', '-k', '5g,5']
                pipeline.append(cat)
                pipeline.append(sort, stdout_path = run.get_single_output_file_for_annotation('features'),
                                stderr_path = run.get_single_output_file_for_annotation('log_stderr'))

                

