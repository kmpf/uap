import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Post_Sawdust(AbstractStep):
    '''
    bla bla
    '''


    def __init__(self, pipeline):
        super(Post_Sawdust, self).__init__(pipeline)

        self.set_cores(2)

        self.add_connection('in/alignments')

        self.add_connection('out/log_stderr')
        self.add_connection('out/alignments')


        self.add_option('library_type', str, choices = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand'], optional=False)
        self.add_option('seq_type', str, choices = ['RNA', 'DNA'],  optional=False)
        self.add_option('read_type', str, choices = ['single', 'paired'], optional=False)
        self.add_option('split_ident', str,  default=' ', optional=False)
        self.add_option('prf', bool,  default=False, optional=True)


        self.require_tool('post_sawdust')
        self.require_tool('samtools')
        self.require_tool('cat')
        self.require_tool('pigz')



    def runs(self, run_ids_connections_files):
        # Compile the list of options


        for run_id in run_ids_connections_files.keys():
            # Check input files
            alignments = run_ids_connections_files[run_id]['in/alignments'][0]
            input_paths = alignments

            with self.declare_run(run_id) as run:
                stdout_path = run.add_output_file(
                    'alignments',
                    '%s-postsawdust.bam' % run_id,
                    [input_paths])


                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-log_stderr.txt' % run_id,
                    [input_paths])



                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
                        samtools_front = [self.get_tool('samtools'), 'view', '-h', alignments]
                        pipe.add_command(samtools_front)

                        post_sawdust = [ self.get_tool('post_sawdust'),
                                         '--library-type', self.get_option('library_type'),
                                         '--seq-type', self.get_option('seq_type'),
                                         '--read-type', self.get_option('read_type'),
                                         '--split-identifier-by', self.get_option('split_ident')]

                        if self.get_option('prf'):
                            post_sawdust.append('--prf')

                        pipe.add_command(post_sawdust, stderr_path=log_stderr)

                        samtools_end = [self.get_tool('samtools'), 'view', '-Shb', '-']

                        pipe.add_command(samtools_end, stdout_path=stdout_path)
