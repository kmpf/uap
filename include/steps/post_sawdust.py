import sys
from abstract_step import *
import process_pool
import yaml


class Post_Sawdust(AbstractStep):

    def __init__(self, pipeline):
        super(Post_Sawdust, self).__init__(pipeline)
        
        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log_stderr')

        self.add_option('library_type', str, choices = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand'])
        self.add_option('seq_type', str, choices = ['RNA', 'DNA'])

        self.require_tool('post_sawdust')
        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cat4m')
        
    def declare_runs(self):

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file., but got this %s" % input_paths)
                run.add_private_info('in-alignment', input_paths[0])
                run.add_output_file('alignments', '%s-fixed.bam' % run_id, input_paths)
                run.add_output_file('log_stderr', '%s-log_stderr.txt' % run_id, input_paths)


                if is_paired_end: 
                    run.add_private_info('read_type', 'paired')
                else:
                    run.add_private_info('read_type', 'single')
                
                    
    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                alignments_path = run.get_private_info('in-alignment')
                cat4m = [self.get_tool('cat4m'), alignments_path]
                samtools_front = [self.get_tool('samtools'), 'view', '-h', '-'] 
                post_sawdust = [self.get_tool('post_sawdust'),
                                '--library-type', self.get_option('library_type'),
                                '--seq-type', self.get_option('seq_type'),
                                '--read-type',run.get_private_info('read_type')]
                samtools = [self.get_tool('samtools'), 'view', '-Shbo',   run.get_single_output_file_for_annotation('alignments'), '-']

                pipeline.append(cat4m)
                pipeline.append(samtools_front)
                pipeline.append(post_sawdust, stderr_path = run.get_single_output_file_for_annotation('log_stderr'))
                pipeline.append(samtools)

