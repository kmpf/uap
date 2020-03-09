import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Piranha(AbstractStep):
    '''
    Piranha is a peak-caller for CLIP- and RIP-Seq data.
    It takes input in BED or BAM format and identifies regions of statistically significant read enrichment.
    Additional covariates may optionally be provided to further inform the peak-calling process.

    This is for the release: Piranha 1.2.0

    http://smithlabresearch.org/software/piranha/
    '''

    def __init__(self, pipeline):
        super(Piranha, self).__init__(pipeline)

        self.set_cores(1)

        # .bed files from mapping tools
        self.add_connection('in/features')

        # .bed files containing piranha output
        # this it what is provided via -o to Piranha
        self.add_connection('out/features')

        self.require_tool('piranha')
        self.require_tool('bamtools')

        # add options
        self.add_option('sort', bool, optional=True, default=False,
                        description="indicates that input is unsorted and "
                        "Piranha should sort it for you"
                        )
        self.add_option('bin_size_reponse', int, optional=True,
                        description="indicates that the response (first "
                        "input file) is raw reads and should be binned "
                        "into bins of this size"
                        )
        self.add_option('p_threshold', float, optional=True,
                        description="significance threshold for sites"
                        )
        self.add_option(
            'dist',
            str,
            optional=True,
            choices=[
                'Poisson',
                'NegativeBinomial',
                'ZeroTruncatedPoisson',
                'ZeroTruncatedNegativeBinomial',
                'PoissonRegression',
                'NegativeBinomialRegression',
                'ZeroTruncatedPoissonRegression',
                'ZeroTruncatedNegativeBinomialRegression'],
            description="Distribution type. Currently supports Poisson,"
            "NegativeBinomial, ZeroTruncatedPoisson, "
            "ZeroTruncatedNegativeBinomial (default with no "
            "covariates), PoissonRegression, NegativeBinomialRegression, "
            "ZeroTruncatedPoissonRegression, "
            "ZeroTruncatedNegativeBinomialRegression "
            "(default with covariates)")

    def runs(self, run_ids_connections_files):
        options = ['sort', 'p_threshold', 'bin_size_reponse', 'dist']

        set_options = [option for option in options if
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
                input_paths = run_ids_connections_files[run_id]['in/features']

                output_file = run.add_output_file(
                    'features', '%s-piranha.out.bed' %
                    run_id, input_paths)

                piranha = [self.get_tool('piranha'), '-o', output_file]
                piranha.extend(option_list)
                piranha.append(input_paths[0])

                with run.new_exec_group() as exec_group:
                    exec_group.add_command(piranha)
