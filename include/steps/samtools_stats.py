import sys
import yaml

from ..abstract_step import *
from .. import process_pool

class SamtoolsStats(AbstractStep):
    '''
    The samtools_stats step can be used to collect read statistics from BAM files using samtools stats.
    '''

    def __init__(self, pipeline):
        super(SamtoolsStats, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/alignments')
        self.add_connection('out/stats')
                
        self.require_tool('cat')
        self.require_tool('samtools')
        self.require_tool('pigz')
        
    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignment file.")

                basename = os.path.basename(input_paths[0]).split('.')[0]
                run.add_output_file('stats', basename + '.bam.stats', input_paths)
                run.add_private_info('in-bam', input_paths[0])

    def execute(self, run_id, run):
        bam_in_path = run.get_private_info('in-bam')
        sam_out_path = run.get_single_output_file_for_annotation('stats')

        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat = [self.get_tool('cat'), bam_in_path]
                samtools = [self.get_tool('samtools'), 'stats']

                pipeline.append(cat)
                pipeline.append(samtools, stdout_path = sam_out_path)

