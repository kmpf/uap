import sys
from abstract_step import AbstractStep

class Preseq(AbstractStep):
    '''
    The preseq package is aimed at predicting the yield of distinct reads from a
    genomic library from an initial sequencing experiment. The estimates can then
    be used to examine the utility of further sequencing, optimize the sequencing
    depth, or to screen multiple libraries to avoid low complexity samples.
    '''

    def __init__(self, pipeline):
        super(Preseq, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/complexity_curve')
        self.add_connection('out/future_yield')
        
        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('preseq')

        # General options
        self.add_option('pe', bool, optional = False, description =
                        'input is paired end read file')
        self.add_option('hist', bool, optional = True, default = False,
                        description = 'input is a text file containing the '
                        'observed histogram')
        self.add_option('vals', bool, optional = True, default = False,
                        description = 'input is a text file containing only '
                        'the observed counts')
        self.add_option('bam', bool, optional = True, default = True,
                        description = 'input is in BAM format')

        # c_curve specific options
        self.add_option('c_curve_step', bool, optional = False, description =
                        'step size in extrapolations (default: 1e+06)')
        self.add_option('c_curve_seg_len', int, optional = True,
                        description = 'maximum segment length when merging '
                        'paired end bam reads (default: 5000)')

        # lc_extrap specific options
        self.add_option('lc_extrap_', bool, optional = True, description =
                        'step size gin extrapolations (default: 1e+06)')
        self.add_option('lc_extrap_', int, optional = True,
                        description = 'maximum segment length when merging '
                        'paired end bam reads (default: 5000)')

        # gc_extrap specific options
        self.add_option('gc_extrap_max_width', int, optional = True,
                        description = 'max fragment length, set equal to read '
                        'length for single end reads')
        self.add_option('gc_extrap_bin_size', int, optional = True,
                        description = 'bin size (default: 10)')
        self.add_option('gc_extrap_extrap', int, optional = True,
                        description = 'maximum extrapolation in base pairs '
                        '(default: 1e+12)')
        self.add_option('gc_extrap_bootstraps', int, optional = True,
                        description = 'number of bootstraps (default: 100)')
        self.add_option('gc_extrap_cval', float, optional = True,
                        description = 'level for confidence intervals (default:'
                        ' 0.95)')
        self.add_option('gc_extrap_terms', int, optional = False,
                        description = 'maximum number of terms')

    def runs(self, run_ids_connections_files):
        general_opts = ['pe', 'hist', 'vals', 'bam']
        c_curve_opts = ['c_curve_step', 'c_curve_seg_len']
        lc_extrap_opts = []
        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]
        option_list = list()
        for option in set_options:
            option_list.append('-%s' % option)
            if not isinstance(self.get_option(option), bool):
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("complexity_curve")
                    run.add_empty_output_connection("future_yield")
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                    with run.new_exec_group() as cc_group:
                        c_curve_out = run.add_output_file(
                            'complexity_curve',
                            '%s_complexity_output.txt' % run_id,
                            input_paths
                        )
                        c_curve = [self.get_tool('preseq'), 'c_curve']
                        c_curve.extend(option_list)
                        c_curve.extend(['-o', c_curve_out, input_paths[0]])
                        cc_group.add_command(c_curve)

                    with run.new_exec_group() as lc_group:
                        lc_extrap_out = run.add_output_file(
                            'future_yield',
                            '%s_future_yield.txt' % run_id,
                            input_paths
                        )
                        lc_extrap = [self.get_tool('preseq'), 'lc_extrap']
                        lc_extrap.extend(option_list)
                        lc_extrap.extend(['-o', c_curve_out, input_paths[0]])
                        lc_group.add_command(lc_extrap)

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
                                        '--processes', '1',
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
