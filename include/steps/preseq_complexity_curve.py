from uaperrors import UAPError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PreseqComplexityCurve(AbstractStep):
    '''
    The preseq package is aimed at predicting the yield of distinct reads from a
    genomic library from an initial sequencing experiment. The estimates can then
    be used to examine the utility of further sequencing, optimize the sequencing
    depth, or to screen multiple libraries to avoid low complexity samples.

    c_curve computes the expected yield of distinct reads for experiments smaller
    than the input experiment in a .bed or .bam file through resampling. The full
    set of parameters can be outputed by simply typing the program name. If
    output.txt is the desired output file name and input.bed is the input .bed
    file, then simply type::

        preseq c_curve -o output.txt input.sort.bed

    Documentation::

        http://smithlabresearch.org/software/preseq/

    '''

    def __init__(self, pipeline):
        super(PreseqComplexityCurve, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/complexity_curve')

        self.require_tool('preseq')

        # c_curve specific options
        self.add_option('step', int, optional = True, description =
                        'step size in extrapolations (default: 1e+06)')
        self.add_option('verbose', bool, optional = True, description =
                        'print more information')
        self.add_option('pe', bool, optional = False, description =
                        'input is paired end read file')
        self.add_option('hist', bool, optional = True, default = False,
                        description = 'input is a text file containing the '
                        'observed histogram')
        self.add_option('vals', bool, optional = True, default = False,
                        description = 'input is a text file containing only '
                        'the observed counts')
        self.add_option('seg_len', int, optional = True, description =
                        'maximum segment length when merging paired end bam '
                        'reads (default: 5000)')

    def runs(self, run_ids_connections_files):
        options = ['step', 'verbose', 'pe', 'hist', 'vals', 'seg_len']

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
                elif len(input_paths) != 1:
                    raise UAPError("Expected exactly one alignments file.")
                elif not is_bam and not is_bed:
                    raise UAPError("Input file %s is niether BAM nor BED." %
                                 input_paths[0])
                else:
                    with run.new_exec_group() as cc_group:
                        c_curve_out = run.add_output_file(
                            'complexity_curve',
                            '%s_complexity_output.txt' % run_id,
                            input_paths
                        )
                        c_curve = [self.get_tool('preseq'), 'c_curve']
                        c_curve.extend(option_list)
                        if is_bam:
                            c_curve.append('-bam')
                        c_curve.extend(['-o', c_curve_out, input_paths[0]])
                        cc_group.add_command(c_curve)
