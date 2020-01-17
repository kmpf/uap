from uaperrors import UAPError
from logging import getLogger
from abstract_step import AbstractStep


logger = getLogger('uap_logger')


class Kallisto(AbstractStep):
    '''
    Kallisto


    '''

    def __init__(self, pipeline):
        super(Kallisto, self).__init__(pipeline)

        # input connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('in/kallisto-index')

        # output connections
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

        self.add_option('index', str, optional=True, default=None,
                        description="Filename for the kallisto index to be \
                        used for quantification")

        self.add_option('bias', bool, optional=True, default=None,
                        description="Perform sequence based bias correction")

        self.add_option('bootstrap-samples', int, optional=True, default=None,
                        description="Number of bootstrap samples (default: 0)")

        self.add_option('seed', int, optional=True, default=None,
                        description="Seed for the bootstrap sampling")

        self.add_option('single', bool, optional=True, default=None,
                        description="Quantify single-end reads")

        self.add_option('single-overhang', bool, optional=True, default=None,
                        description="""Include reads where unobserved
                        rest of fragment is predicted to lie
                        outside a transcript""")

        self.add_option('fr-stranded', bool, optional=True, default=None,
                        description="Strand specific reads, first read forward")

        self.add_option('rf-stranded', bool, optional=True, default=None,
                        description="Strand specific reads, first read reverse")

        self.add_option('fragment-length', int, optional=True, default=None,
                        description="Estimated average fragment length")

        self.add_option('sd', int, optional=True, default=None,
                        description="Estimated standard deviation of fragment length")
        # skipping gtf, pseudobam, genomebam, chromosomes

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
                input_fileset = []
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                input_fileset.append(r1)

                r2 = None
                if 'in/second_read' in run_ids_connections_files[run_id]:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]
                    input_fileset.append(r2)

                kallisto_eg = run.new_exec_group()

                kallisto = [self.get_tool('kallisto'), 'quant']
                kallisto.extend(['-t', str(self.get_option('cores'))])

                d_files = input_fileset[:]
                if self.is_option_set_in_config('index'):
                    kallisto.extend(['--index', self.get_option('index')])
                else:
                    if connect_index_path:
                        kallisto.extend(['--index', connect_index_path])
                        d_files.append( connect_index_path)
                    else:
                        raise UAPError(
                        "%s no kallisto index give via config or connection" % run_id)

                flags = ['fr-stranded', 'rf-stranded',
                         'bias', 'single-overhang', 'single']

                for flag in flags:
                    if self.is_option_set_in_config(flag) and self.get_option(flag):
                        kallisto.append('--' + flag)

                param_flags = ['bootstrap-samples', 'seed',
                               'fragment-length', 'sd']

                for param_flag in param_flags:
                    if self.is_option_set_in_config(param_flag):
                        kallisto.extend(['--' + param_flag,
                                         str(self.get_option(param_flag))])

                out_path = run.get_output_directory_du_jour_placeholder()
                kallisto.extend(['-o', out_path])

                kallisto.extend(input_fileset)

                stderr_file = "%s-kallisto-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, d_files)
                stdout_file = "%s-kallisto-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, d_files)

                h5_file = "abundance.h5"
                run.add_output_file("abundance.h5", h5_file, d_files)
                tsv_file = "abundance.tsv"
                run.add_output_file("abundance.tsv", tsv_file, d_files)
                run_info = "run_info.json"
                run.add_output_file("run_info.json", run_info, d_files)

                kallisto_eg.add_command(kallisto, stdout_path=log_stdout,
                                        stderr_path=log_stderr)
