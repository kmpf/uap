import sys
from abstract_step import *
import process_pool
import yaml


class BamToSam(AbstractStep):

    def __init__(self, pipeline):
        super(BamToSam, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('pigz')

        
        self.add_option('genome', str, optional=False)

    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")

                basename = os.path.basename(input_paths[0]).split('.')[0]
                run.add_output_file('alignments', basename + '.sam.gz', input_paths)
                run.add_private_info('in-bam', input_paths[0])
                

    def execute(self, run_id, run):
        bam_in_path = run.get_private_info('in-bam')
        sam_out_path = run.get_single_output_file_for_annotation('alignments')
        

        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), bam_in_path]
                samtools = [self.get_tool('samtools'), 'view', '-h', '-']
                pigz2 = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']

                pipeline.append(cat4m)
                pipeline.append(samtools)
                pipeline.append(pigz2, stdout_path = sam_out_path)

