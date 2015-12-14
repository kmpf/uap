import sys
import os
from abstract_step import AbstractStep

class PreseqFutureYield(AbstractStep):
    '''
    The preseq package is aimed at predicting the yield of distinct reads from a
    genomic library from an initial sequencing experiment. The estimates can then
    be used to examine the utility of further sequencing, optimize the sequencing
    depth, or to screen multiple libraries to avoid low complexity samples.
    '''

    def __init__(self, pipeline):
        super(PreseqFutureYield, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/future_yield')
        
        self.require_tool('preseq')

        # lc_extrap specific options
        self.add_option('extrap', int, optional = True, description =
                        'maximum extrapolation (default: 1e+10)')
        self.add_option('step', int, optional = True, description =
                        'step size in extrapolations (default: 1e+06)')
        self.add_option('bootstraps', int, optional = True, description =
                        'number of bootstraps (default: 100)')
        self.add_option('cval', float, optional = True, description =
                        'level for confidence intervals (default: 0.95)')
        self.add_option('dupl_level', float, optional = True, description =
                        'fraction of duplicate to predict (default: 0.5)')
        self.add_option('terms', int, optional = True, description =
                        'maximum number of terms')
        self.add_option('seg_len', int, optional = True, description =
                        'maximum segment length when merging paired end bam '
                        'reads (default: 5000)')
        self.add_option('pe', bool, optional = False, description =
                        'input is paired end read file')
        self.add_option('hist', bool, optional = True, default = False,
                        description = 'input is a text file containing the '
                        'observed histogram')
        self.add_option('vals', bool, optional = True, default = False,
                        description = 'input is a text file containing only '
                        'the observed counts')
        self.add_option('quick', bool, optional = True, description =
                        'quick mode, estimate yield without bootstrapping for '
                        'confidence intervals')

    def runs(self, run_ids_connections_files):
        options = ['extrap', 'step', 'bootstraps', 'cval', 'dupl_level',
                        'terms', 'seg_len', 'pe', 'hist', 'vals', 'quick']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool) and \
               self.get_option(option):
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
                    run.add_empty_output_connection("future_yield")
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                elif not is_bam and not is_bed:
                    raise StandardError("Input file %s is niether BAM nor BED." %
                                        input_paths[0])
                else:
                    with run.new_exec_group() as lc_group:
                        lc_extrap_out = run.add_output_file(
                            'future_yield',
                            '%s_future_yield.txt' % run_id,
                            input_paths
                        )
                        lc_extrap = [self.get_tool('preseq'), 'lc_extrap']
                        lc_extrap.extend(option_list)
                        if is_bam:
                            lc_extrap.append('-bam')
                        lc_extrap.extend(['-o', lc_extrap_out, input_paths[0]])
                        lc_group.add_command(lc_extrap)

