import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Reformatcigar(AbstractStep):
    '''
    The segemehl release from April 2017 (segemehl/export.13) uses the 
    expanded version of Cigar strings: It clearly separates between match (=)
    and mismatch (X) while in older Cigar string definitions both are 
    encoded with an M character.
    However, subsequent tools like htseq-count are not able to handle these 
    newer versions of Cigar strings.
    
    This step transforms the new Cigar string encoding back to the old one.

    This step depends on the segemehl_2017 step. You still may need to 
    add the subsequenct step s2c before calling cufflinks and/or htseq-count.
    '''

    def __init__(self, pipeline):
        super(Reformatcigar, self).__init__(pipeline)
        
        # connections - indentifier for in/output
        #             - expects list, maybe empty or 'none', 
        #               e.g. if only first_read info available
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log')
        
        # external tools
        self.require_tool('pigz')
        self.require_tool('cat')
        # internal tools
        self.require_tool('segemehl_2017_to_reformatCigar')

        # step options
        self.add_option('threads', int, optional=True,
                        description='Number of threads 2B started. (Default: 1). '
                        'Beware that this is only for (un-)compressing, '
                        'the reformating is using a single CPU only.')

    def runs(self, run_ids_connections_files):

        # set # of cores for cluster, it is ignored if run locally
        if self.get_option('threads'):
            self.set_cores(self.get_option('threads'))
        else:
            self.set_cores(1)

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:

                # get input file name for each run
                input_paths = run_ids_connections_files[run_id]['in/alignments']

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file."
                                        ", but got this %s" % input_paths)

                alignments_path = input_paths[0]

                cat = [self.get_tool('cat'), alignments_path]
                pigzD = [self.get_tool('pigz'), '--decompress',
                         '--processes', str(self.get_cores()),
                         '--stdout']
                reformatcigar = [self.get_tool('segemehl_2017_to_reformatCigar'),
                                 '--in-file', '/dev/stdin']
                pigzC = [self.get_tool('pigz'), 
                         '--processes', str(self.get_cores()),
                         '--stdout']

                out_file = run.add_output_file('alignments',
                                               '%s-reformatCigar.sam.gz' % run_id,
                                               input_paths)
                log_file = run.add_output_file('log', 
                                               '%s-reformatCigar.err.log' % run_id,
                                               input_paths)

                with run.new_exec_group() as exec_group:
                    
                    with exec_group.add_pipeline() as reformatCigar_pipe:

                        reformatCigar_pipe.add_command(cat)
                        reformatCigar_pipe.add_command(pigzD)
                        reformatCigar_pipe.add_command(reformatcigar)
                        reformatCigar_pipe.add_command(pigzC,
                                                       stdout_path = out_file,
                                                       stderr_path = log_file)
                        
