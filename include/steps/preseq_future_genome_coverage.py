import sys
from logging import getLogger
from abstract_step import AbstractStep
import os

logger=getLogger('uap_logger')

class PreseqFutureGenomeCoverage(AbstractStep):
    '''
    The preseq package is aimed at predicting the yield of distinct reads from a
    genomic library from an initial sequencing experiment. The estimates can then
    be used to examine the utility of further sequencing, optimize the sequencing
    depth, or to screen multiple libraries to avoid low complexity samples.

    gc_extrap computes the expected genomic coverage for deeper sequencing for
    single cell sequencing experiments. The input should be a mr or bed file.
    The tool bam2mr is provided to convert sorted bam or sam files to mapped
    read format.
    '''

    def __init__(self, pipeline):
        super(PreseqFutureGenomeCoverage, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/future_genome_coverage')

        self.require_tool('preseq')

        # gc_extrap specific options

        self.add_option('max_width', int, optional = True, description =
                        'max fragment length, set equal to read length for '
                        'single end reads')
        self.add_option('bin_size', int, optional = True,
                        description = 'bin size (default: 10)')
        self.add_option('extrap', int, optional = True, description =
                        'maximum extrapolation in base pairs (default: 1e+12)')
        self.add_option('step', int, optional = True, description =
                        'step size in bases between extrapolations (default: '
                        '1e+08)')
        self.add_option('bootstraps', int, optional = True, description =
                        'number of bootstraps (default: 100)')
        self.add_option('cval', float, optional = True, description =
                        'level for confidence intervals (default: 0.95)')
        self.add_option('terms', int, optional = True, description =
                        'maximum number of terms')
        self.add_option('quick', bool, optional = True, description =
                        'quick mode: run gc_extrap without bootstrapping for '
                        'confidence intervals')

    def runs(self, run_ids_connections_files):
        options = ['max_width', 'bin_size', 'extrap', 'step', 'bootstraps',
                   'cval', 'terms', 'quick']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                is_bam = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.bam'] else False
                is_bed = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.bed'] else False

                if input_paths == [None]:
                    run.add_empty_output_connection("complexity_curve")
                    run.add_empty_output_connection("future_yield")
                elif len(input_paths) != 1:
                    logger.error("Expected exactly one alignments file.")
                    sys.exit(1)
                elif not is_bam and not is_bed:
                    logger.error("Input file %s is niether BAM nor BED." %
                                 input_paths[0])
                    sys.exit(1)
                else:
                    with run.new_exec_group() as gc_group:
                        gc_extrap_out = run.add_output_file(
                            'future_genome_coverage',
                            '%s_future_genome_coverage.txt' % run_id,
                            input_paths
                        )
                        gc_extrap = [self.get_tool('preseq'), 'gc_extrap']
                        gc_extrap.extend(option_list)
                        if is_bed:
                            gc_extrap.append('-bed')
                        gc_extrap.extend(['-o', gc_extrap_out, input_paths[0]])
                        gc_group.add_command(gc_extrap)
