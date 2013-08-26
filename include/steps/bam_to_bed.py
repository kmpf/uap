import sys
from abstract_step import *
import process_pool
import re
import yaml


class BamToBed(AbstractStep):

    def __init__(self, pipeline):
        super(BamToBedgraph, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('bedtools')
        self.require_tool('mate_pair_strand_switch')

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignment file. Got %s"
                                        % input_paths)
                if  input_paths[0][-4:] != '.bam':
                    raise StandardError("%s file suffix is not '.bam'. " +
                                        "Please provide a BAM file" % input_paths[0])
                
                bed_file = os.path.basename(input_paths[0])[:-4] + '.bed'
                run.add_output_file('alignments', bed_file, input_paths)
                run.add_private_info('out-bed', bed_file)
                run.add_private_info('in-bam', input_paths[0])
                
    def execute(self, run_id, run):
        bam_path = run.get_private_info('in-bam')
        bed_file = run.get_private_info('out-bed')
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m_in = [self.get_tool('cat4m'), bam_path]
                samtools = [self.get_tool('samtools'), '-bf', '0x2', '-']
                bedtools = [self.get_tool('bedtools'), 'bamtobed', '-i', 'stdin']
                strand_switch = [self.get_tool('mate_pair_strand_switch')]
                
                pipeline.append(cat4m)
                pipeline.append(samtools)
                pipeline.append(bedtools)
                pipeline.append(strand_switch, stdout_path = bed_file)
