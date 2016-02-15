import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class BamToBedgraph(AbstractStep):

    def __init__(self, pipeline):
        super(BamToBedgraph, self).__init__(pipeline)
        
        self.set_cores(8)
        
        self.add_connection('in/alignments',
                            constraints = {'min_files_per_run': 1,
                                           'max_files_per_run': 1})
        self.add_connection('out/alignments')
        
        self.require_tool('dd')
        self.require_tool('pigz')
        self.require_tool('mkfifo')
        self.require_tool('bedtools')
        self.require_tool('bedToBigBed')
        self.require_tool('bedGraphToBigWig')

        # General Options
        self.add_option('output-format', str,
                        choices = ['bed', 'bigBed', 'bedGraph', 'bigWig'],
                        default = 'bigWig',
                        optional = False)
        self.add_option('chromosome-sizes', str, optional = False)
        # Options for bedtools bamtobed
        self.add_option('bedtools-bamtobed-tag', str, optional = True)
        self.add_option('bedtools-bamtobed-color', str, optional = True)

        # Options for bedtools genomecov (that make sense for BAM to BG)
        self.add_option('bedtools-genomecov-report-zero-coverage',
                        bool, optional = False)
        self.add_option('bedtools-genomecov-max', int, optional = True)
        self.add_option('bedtools-genomecov-split', bool, default = True)
        self.add_option('bedtools-genomecov-strand', str, choices = ['+', '-'],
                        optional = True)
        self.add_option('bedtools-genomecov-scale', float, optional = True)
        self.add_option('bedtools-genomecov-5', bool,
                        default = False, optional = True)
        self.add_option('bedtools-genomecov-3', bool,
                        default = False, optional = True)


        self.add_option('trackline', dict, optional = True)
        self.add_option('trackopts', dict, optional = True)

    def runs(self, run_ids_connections_files):
        def compile_option_list(prefix, options):
            options = ['%s%s' % (prefix, x) for x in options]
            set_options = [option for option in options if \
                           self.is_option_set_in_config(option)]

            option_list = list()
            for option in set_options:
                if isinstance(self.get_option(option), bool):
                    if self.get_option(option):
                        option_list.append('-%s' % option.replace(prefix, ''))
                else:
                    option_list.append('-%s' % option.replace(prefix, ''))
                    option_list.append(str(self.get_option(option)))
            return option_list

        # Compile the list of options for bedtools genomecov
        bedtools_genomecov_options = compile_option_list(
            'bedtools-genomecov-', ['max', 'split', '5', '3', 'strand', 'scale']
        )
        bedtools_bamtobed_options = compile_option_list(
            'bedtools-bamtobed-', ['tag', 'color']
        )

        # Check if chromosome sizes points to a real file
        if not os.path.isfile(self.get_option('chromosome-sizes')):
            logger.error("Value for option 'chromosome-sizes' is not a "
                         "file: %s" % self.get_option('chromosome-sizes'))
            sys.exit(1)
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                is_gzipped = False
                is_bam = False
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    logger.error("Expected exactly one alignments file.")
                    sys.exit(1)
                else:
                    root, ext = os.path.splitext(os.path.basename(input_paths[0]))
                    if ext in ['.gz', '.gzip']:
                        is_gzipped = True
                    elif ext in ['.bam']:
                        is_bam = True
                    if is_gzipped:
                        root, ext = os.path.splitext(root)
                        if ext in ['.bam']:
                            is_bam = True

                    if not is_bam:
                        logger.error("The file %s does not appear to be any "
                                     "of bam.gz, bam.gzip, or bam"
                                     % input_paths[0]
                        )
                        sys.exit(1)
                    with run.new_exec_group() as exec_group:
                        # Create FIFO for use with bedToBigBed, bedGraphToBigWig
                        big_fifo = run.add_temporary_file(
                            'big_fifo', designation = 'ouput')
                        if self.get_option('output-format') in \
                           ['bigBed', 'bigWig']:
                            mkfifo = [self.get_tool('mkfifo'), big_fifo]
                            exec_group.add_command(mkfifo)

                        with exec_group.add_pipeline() as pipe:
                            # 1. command: Read file in 4MB chunks
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

                            output_file = str()
                            
                            # BAM -> BED
                            # (necessary for bed, bigBed)
                            if self.get_option('output-format') in \
                               ['bed', 'bigBed']:
                                bam_to_bed = [
                                    self.get_tool('bedtools'), 'bamtobed']
                                bam_to_bed.extend(bedtools_bamtobed_options)
                                bam_to_bed.extend(['-i', 'stdin'])

                                pipe.add_command(bam_to_bed)
                                # Set output file name here for dd_out
                                output_file = '%s.bed' % root

                            # BAM -> BedGraph
                            # (necessary for bedGraph, bigWig)
                            if self.get_option('output-format') in \
                               ['bedGraph', 'bigWig']:
                                # 2. command: Convert BAM to BedGraph
                                genomecov = [
                                    self.get_tool('bedtools'), 'genomecov']
                                genomecov.extend( bedtools_genomecov_options)
                                genomecov.extend( ['-ibam', 'stdin'] )

                                if self.get_option('bedtools-genomecov-report-zero-coverage'):
                                    genomecov.append('-bga')
                                else:
                                    genomecov.append('-bg')

                                pipe.add_command(genomecov)
                                # Set output file name here for dd_out
                                output_file = '%s.bg' % root
                                    
                            # Write BED or BedGraph output to file
                            if self.get_option('output-format') in \
                               ['bed', 'bedGraph']:
                                dd_out = [self.get_tool('dd'), 'obs=4M']
                                pipe.add_command(
                                    dd_out,
                                    stdout_path = run.add_output_file(
                                        'alignments',
                                        output_file,
                                        input_paths)
                                )

                            # BED -> BigBed
                            if self.get_option('output-format') in ['bigBed']:
                                bed_to_bigbed = [
                                    self.get_tool('bedToBigBed'),
                                    'stdin',
                                    self.get_option('chromosome-sizes'),
                                    big_fifo
                                ]
                                output_file = '%s.bb' % root
                                pipe.add_command(bed_to_bigbed)
                            
                            # BedGraph -> BigWig
                            if self.get_option('output-format') in ['bigWig']:
                                bedgraph_to_bigwig = [
                                    self.get_tool('bedGraphToBigWig'),
                                    'stdin',
                                    self.get_option('chromosome-sizes'),
                                    big_fifo
                                ]
                                output_file = '%s.bw' % root
                                pipe.add_command(bedgraph_to_bigwig)

                            if self.get_option('output-format') in \
                               ['bedGraph', 'bigWig']:
                                dd_out = [self.get_tool('dd'),
                                          'bs=4M',
                                          'if=%s' % big_fifo]

                                pipe.add_command(
                                    dd_out,
                                    stdout_path = run.add_output_file(
                                        'alignments',
                                        output_file,
                                        input_paths)
                                )
