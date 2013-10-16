import sys
from abstract_step import *
import process_pool
import yaml


class SamToBam(AbstractStep):

    def __init__(self, pipeline):
        super(SamToBam, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/indices')
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('pigz')

        self.add_option('sort_by_name', bool, default= False)
        self.add_option('genome', str, optional=False)

    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                
                basename = os.path.basename(input_paths[0]).split('.')

                if 'sam' in basename:
                    sam_index = basename.index('sam')
                    basename = basename[:sam_index]
                elif 'bam' in basename:
                    bam_index  = basename.index('bam')
                    basename = basename[:bam_index]
                else:
                    raise StandardError("File %s is neither a BAM nor SAM file" % (input_paths[0]))
                basename = '.'.join(basename)
                run.add_output_file('alignments', basename + '.bam', input_paths)
                run.add_output_file('indices', basename + '.bam.bai', input_paths)

                run.add_private_info('in-sam', input_paths[0])

    def execute(self, run_id, run):
        sam_path = run.get_private_info('in-sam')
        sorted_bam_path = run.get_single_output_file_for_annotation('alignments')
        sorted_bai_path = run.get_single_output_file_for_annotation('indices')
        unsorted_bam_path = self.get_temporay_path('sam_to_bam_unsorted', 'output')
        
        use_unsorted_bam_input = unsorted_bam_path
        
        # samtools view

        if sam_path[-7:] == '.sam.gz':
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), sam_path]
                    pigz1 = [self.get_tool('pigz'), '--processes', '2', '-d', '-c']
                    samtools = [self.get_tool('samtools'), 'view', '-Sbt', self.get_option('genome'), '-']
                    
                    pipeline.append(cat4m)
                    pipeline.append(pigz1)
                    pipeline.append(samtools, stdout_path = unsorted_bam_path)
        else:
            # it must be a BAM file already
            use_unsorted_bam_input = sam_path
            
        # samtools sort

        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), use_unsorted_bam_input]
                samtools = [self.get_tool('samtools'), 'sort']
                if self.get_option('sort_by_name'):
                    samtools.append('-n')
                samtools.extend(['-', sorted_bam_path[:-4]])
                
                pipeline.append(cat4m)
                pipeline.append(samtools, hints = {'writes': [sorted_bam_path]})

        # samtools index
        
        with process_pool.ProcessPool(self) as pool:
            pool.launch([self.get_tool('samtools'), 'index', sorted_bam_path, '/dev/stdout'],
                stdout_path = sorted_bai_path, hints = {'reads': [sorted_bam_path]})

        if os.path.exists(unsorted_bam_path):
            os.unlink(unsorted_bam_path)
