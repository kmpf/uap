import sys
from abstract_step import *
import process_pool
import yaml


class SamToSortedBam(AbstractStep):

    def __init__(self, pipeline):
        super(SamToSortedBam, self).__init__(pipeline)
        
        self.set_cores(12)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('pigz')

        self.add_option('sort-by-name', bool, default = False)
        self.add_option('genome-faidx', str, optional = False)


    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                    with run.new_exec_group() as exec_group:

                        with exec_group.add_pipeline() as pipe:
                            # 1 command: Read file in 4MB chunks
                            dd_in = [self.get_tool('dd'),
                                     'ibs=4M',
                                     'if=%s' % input_paths[0]]
                            pipe.add_command(dd_in)

                            if is_gzipped:
                                # 1.1 command: Uncompress file to fifo
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--processes', '1'
                                        '--stdout']
                                pipe.add_command(pigz)

                            # 2 command: Convert sam to bam
                            samtools_view = [
                                self.get_tool('samtools'), 'view',
                                '-S', '-b', '-t',
                                self.get_option('genome-faidx'), '-',
                                '-@', '2'
                            ]
                            pipe.add_command(samtools_view)

                            # 3 command: Sort BAM input
                            samtools_sort = [
                                self.get_tool('samtools'), 'sort',
                                '-O', 'bam'
                            ]
                            if self.get_option('sort-by-name'):
                                samtools.append('-n')
                            samtools_sort.extend(
                                ['-T', run_id, 
                                 '-',
                                 '-@', '6']
                            )
                            pipe.add_command(samtools_sort)

                            # 4 command:
                            dd_out = [self.get_tool('dd'), 'obs=4M']
                            pipe.add_command(
                                dd_out,
                                stdout_path = run.add_output_file(
                                    'alignments',
                                    run_id + '.sorted.bam',
                                    input_paths)
                            )


#                            # 4 command: Index sorted BAM
#                            samtools_index = 
#            pool.launch([self.get_tool('samtools'), 'index', sorted_bam_path, '/dev/stdout'],
#                stdout_path = sorted_bai_path, hints = {'reads': [sorted_bam_path]})
#
#
#                basename = os.path.basename(input_paths[0]).split('.')
#
#                if 'sam' in basename:
#                    sam_index = basename.index('sam')
#                    basename = basename[:sam_index]
#                elif 'bam' in basename:
#                    bam_index  = basename.index('bam')
#                    basename = basename[:bam_index]
#                else:
#                    raise StandardError("File %s is neither a BAM nor SAM file" % (input_paths[0]))
#                basename = '.'.join(basename)
#                run.add_output_file('alignments', basename + '.bam', input_paths)
#                run.add_output_file('indices', basename + '.bam.bai', input_paths)
#
#                run.add_private_info('in-sam', input_paths[0])
#                run.new_exec_group()
#
#    def execute(self, run_id, run):
#        sam_path = run.get_private_info('in-sam')
#        sorted_bam_path = run.get_single_output_file_for_annotation('alignments')
#        sorted_bai_path = run.get_single_output_file_for_annotation('indices')
#        unsorted_bam_path = self.get_temporary_path('sam_to_bam_unsorted', 'output')
#        
#        use_unsorted_bam_input = unsorted_bam_path
#        
#        # samtools view
#
#        if sam_path[-7:] == '.sam.gz':
#            with process_pool.ProcessPool(self) as pool:
#                with pool.Pipeline(pool) as pipeline:
#                    cat = [self.get_tool('cat'), sam_path]
#                    pigz1 = [self.get_tool('pigz'), '--processes', '2', '-d', '-c']
#                    
#                    
#                    pipeline.append(cat)
#                    pipeline.append(pigz1)
#                    pipeline.append(samtools, stdout_path = unsorted_bam_path)
#        else:
#            # it must be a BAM file already
#            use_unsorted_bam_input = sam_path
#            
#        # samtools sort
#
#        with process_pool.ProcessPool(self) as pool:
#            with pool.Pipeline(pool) as pipeline:
#                cat = [self.get_tool('cat'), use_unsorted_bam_input]
#                samtools = [self.get_tool('samtools'), 'sort']
#                if self.get_option('sort_by_name'):
#                    samtools.append('-n')
#                samtools.extend(['-', sorted_bam_path[:-4]])
#                
#                pipeline.append(cat)
#                pipeline.append(samtools, hints = {'writes': [sorted_bam_path]})
#
#        # samtools index
#        
#        with process_pool.ProcessPool(self) as pool:
#            pool.launch([self.get_tool('samtools'), 'index', sorted_bam_path, '/dev/stdout'],
#                stdout_path = sorted_bai_path, hints = {'reads': [sorted_bam_path]})
