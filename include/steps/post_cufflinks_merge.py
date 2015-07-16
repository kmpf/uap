import sys
from abstract_step import *
import process_pool
import yaml


class Post_Cufflinks_Merge(AbstractStep):

    def __init__(self, pipeline):
        super(Post_Cufflinks_Merge, self).__init__(pipeline)
        
        self.set_cores(2)
        
        self.add_connection('in/features')
        self.add_connection('out/features')
        self.add_connection('out/log_stderr')

        self.add_option('remove_gencode'     , bool,   default=False )
        self.add_option('remove_unstranded'  , bool,   default=False)
        self.add_option('remove_by_gene_name', bool, default=False)
        self.add_option('string', str, optional=True )
        self.add_option('filter_by_class', bool, default=False)
        self.add_option('class_list', str, optional=True)
        self.add_option('filter_by_class_and_gene_name', bool, default=False)

        self.require_tool('post_cufflinks_merge')
        self.require_tool('cat4m')
        
    def declare_runs(self):

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/features'):
            with self.declare_run(run_id) as run:
                run.new_exec_group()
                run.add_private_info('in-features', input_paths[0])
                run.add_output_file('features', '%s-novel.gtf' % run_id, input_paths)
                run.add_output_file('log_stderr', '%s-log_stderr.txt' % run_id, input_paths)


                    
    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                features_path = run.get_private_info('in-features')
                cat4m = [self.get_tool('cat4m'), features_path]

                post_cufflinks_merge = [self.get_tool('post_cufflinks_merge')]


                if self.get_option('remove_gencode'):
                    post_cufflinks_merge.extend(['--remove-gencode'])
                if self.get_option('remove_unstranded'):
                    post_cufflinks_merge.extend(['--remove-unstranded'])
                if self.get_option('remove_by_gene_name'):
                    post_cufflinks_merge.extend(['--remove-by-gene-name'])
                if self.get_option('filter_by_class'):
                    post_cufflinks_merge.extend(['--filter-by-class'])
                if self.get_option('filter_by_class_and_gene_name'):
                    post_cufflinks_merge.extend(['--filter-by-class-and-gene-name'])
                if self.get_option('string'):
                    post_cufflinks_merge.extend(['--string', self.get_option('string')])
                if self.get_option('class_list'):
                    post_cufflinks_merge.extend(['--class-list', self.get_option('class_list')])

                pipeline.append(cat4m)
                pipeline.append(post_cufflinks_merge, stdout_path = run.get_single_output_file_for_annotation('features'),
                                stderr_path = run.get_single_output_file_for_annotation('log_stderr'))
                

