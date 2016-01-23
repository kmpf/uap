import sys
from abstract_step import *

class SamtoolsStats(AbstractStep):
    '''
    The samtools_stats step can be used to collect read statistics from BAM
    files using samtools stats.
    '''

    def __init__(self, pipeline):
        super(SamtoolsStats, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/alignments')
        self.add_connection('out/stats')
                
        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('pigz')

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():
            # Get input alignments
            input_paths = run_ids_connections_files[run_id]\
                          ['in/alignments']
            with self.declare_run(run_id) as run:
                for input_path in input_paths:
                    basename = od.path.splitext(
                        os.path.basename(input_path))[0]

                    with run.new_exec_group().add_pipeline() as pipe:
                        # Read input alignments
                        dd = [self.get_tool('dd'),
                              'ibs=4M',
                              'if=%s' % input_path]
                        pipe.add_command(dd)
                        # Assemble samtools stats command
                        samtools = [self.get_tool('samtools'), 'stats']
                        pipe.add_command(samtools,
                                         stdout_path = run.add_output_file(
                                             'stats',
                                             basename + '.bam.stats',
                                             input_path)
                                     )
