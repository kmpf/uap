from uaperrors import StepError
import sys
from logging import getLogger
from abstract_step import AbstractStep
import os

logger = getLogger('uap_logger')


class Kallisto(AbstractStep):
    '''
    Kallisto with fusion detection option


    '''

    def __init__(self, pipeline):
        super(Kallisto, self).__init__(pipeline)

        # input connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('in/kallisto-index')

        # output connections
        self.add_connection('out/fusion.txt')
        self.add_connection('out/abundance.h5')
        self.add_connection('out/abundance.tsv')
        self.add_connection('out/run_info.json')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        # required tools
        self.require_tool('kallisto')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                        engine and threads ie")

        self.add_option('index', str, optional=False, default=None,
                        description="Filename for the kallisto index to be \
                        used for quantification")

        self.add_option('bias', bool, optional=True, default=None,
                        description="Perform sequence based bias correction")

        self.add_option('bootstrap-samples', int, optional=True, default=None,
                        description="Number of bootstrap samples (default: 0)")

        self.add_option('seed', int, optional=True, default=None,
                        description="Seed for the bootstrap sampling")

        self.add_option('single', int, optional=True, default=None,
                        description="Quantify single-end reads")

        self.add_option('single-overhang', int, optional=True, default=None,
                        description="""Include reads where unobserved
                        rest of fragment is predicted to lie
                        outside a transcript""")

        self.add_option(
            'fr-stranded',
            bool,
            optional=True,
            default=None,
            description="Strand specific reads, first read forward")

        self.add_option(
            'rf-stranded',
            bool,
            optional=True,
            default=None,
            description="Strand specific reads, first read reverse")

        self.add_option('fragment-length', int, optional=True, default=None,
                        description="Estimated average fragment length")

        self.add_option(
            'sd',
            int,
            optional=True,
            default=None,
            description="Estimated standard deviation of fragment length")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        connect_index_path = None
        for run_id in run_ids_connections_files.keys():
            try:
                connect_index_path = run_ids_connections_files[run_id]['in/kallisto-index'][0]
            except KeyError:
                pass

        for run_id in run_ids_connections_files.keys():
            if 'in/kallisto-index' in run_ids_connections_files[run_id]:
                continue

            with self.declare_run(run_id) as run:
                input_path = []
                # Get lst of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read'][0]
                sr_input = run_ids_connections_files[run_id]['in/second_read'][0]

                # Do we have paired end data and is it exactly one ?
                input_paths = [fr_input]

                if sr_input is None:
                    raise StepError(self, "Not paired end")
                else:
                    input_paths.append(sr_input)

                # create logfiles
                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-chimpipe-log_stderr.txt' % run_id,
                    input_paths)

                log_stdout = run.add_output_file(
                    'log_stdout',
                    '%s-chimpipe-log_stdout.txt' % run_id,
                    input_paths)

                run.add_output_file(
                    'fusion.txt', 'fusion.txt', input_paths)

                run.add_output_file(
                    'abundance.h5', 'abundance.h5', input_paths)

                run.add_output_file(
                    'abundance.tsv', 'abundance.tsv', input_paths)

                run.add_output_file(
                    'run_info.json', 'run_info.json', input_paths)

                # Assemble kallisto command
                with run.new_exec_group() as exec_group:

                    kallisto = [self.get_tool('kallisto'),
                                'quant',
                                '-t', str(self.get_option('cores')),
                                '--fusion']
                    if self.is_option_set_in_config('index'):
                        kallisto.extend(
                            ['-i', os.path.abspath(self.get_option('index'))])
                    else:
                        if connect_index_path:
                            kallisto.extend(['-i', connect_index_path])
                        else:
                            raise StepError(
                                self,
                                "%s no kallisto index give via config or connection" %
                                run_id)

                    optns = ['fr-stranded', 'rf-stranded',
                             'bias', 'single-overhang', 'single']

                    for optn in optns:
                        if self.is_option_set_in_config(
                                optn) and self.get_option(optn):
                            kallisto.append('--' + optn)

                    param_optns = ['bootstrap-samples', 'seed',
                                   'fragment-length', 'sd']

                    for param_optn in param_optns:
                        if self.is_option_set_in_config(param_optn):
                            kallisto.extend(['--' + param_optn,
                                             str(self.get_option(param_optn))])

                    kallisto.extend(['-o', '.',
                                     fr_input,
                                     sr_input])

                    exec_group.add_command(kallisto,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)
