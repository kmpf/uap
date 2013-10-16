import sys
from abstract_step import *
import process_pool
import re
import yaml


class BamToBed(AbstractStep):
    '''
    This step converts BAM files into BED files. It executes different command line
    pipelines for single-end and paired-end data.

    Single-end read data is converted via::
    
        cat4m <bam-file> | bedtools bamtobed -i stdin | sort -k1,1 -k2,2n > <bed-file>

    Paired-end read data is converted via::

        cat4m <bam-file> | samtools view -bf 0x2 - | bedtools bamtobed -i stdin |
        mate_pair_strand_switch | sort -k1,1 -k2,2n > <bed-file>

    '''

    def __init__(self, pipeline):
        super(BamToBed, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('bedtools')
        self.require_tool('mate_pair_strand_switch')
        self.require_tool('sort')

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
                run.add_private_info('paired_end', is_paired_end)

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignment file. Got %s"
                                        % input_paths)
                if  input_paths[0][-4:] != '.bam':
                    raise StandardError("%s file suffix is not '.bam'. "
                                        "Please provide a BAM file" % input_paths[0])
                
                bed_file = os.path.basename(input_paths[0])[:-4] + '.bed'
                run.add_output_file('alignments', bed_file, input_paths)
#                print('alignments: %s , %s'  % (bed_file, input_paths))
                run.add_private_info('out-bed', bed_file)
#                print('out-bed: %s' % bed_file)
                run.add_private_info('in-bam', input_paths[0])
#                print('in-bam: %s' % input_paths[0])
#                exit(1)
                
    def execute(self, run_id, run):
        is_paired_end = run.get_private_info('paired_end')
        bam_path = run.get_private_info('in-bam')
        bed_file = run.get_private_info('out-bed')
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), bam_path]
                samtools = [self.get_tool('samtools'), 'view', '-bf', '0x2', '-']
                bedtools = [self.get_tool('bedtools'), 'bamtobed', '-i', 'stdin']
                strand_switch = [self.get_tool('mate_pair_strand_switch')]
                sort = [self.get_tool('sort'), '-k1,1', '-k2,2n']
                
                pipeline.append(cat4m)
                if is_paired_end:
                    pipeline.append(samtools)
                pipeline.append(bedtools)
                if is_paired_end:
                    pipeline.append(strand_switch)
                pipeline.append(sort, stdout_path = run.get_single_output_file_for_annotation('alignments'))
