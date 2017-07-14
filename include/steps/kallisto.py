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

        self.add_option('i', str, optional=False, default=None,
                        description="Filename for the kallisto index to be \
                        used for quantification")

        #self.add_option('o', str, optional=False, default=None,
        #                description="Directory to write output to")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                print(run_id)
                input_fileset = []
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                input_fileset.append(r1)

                r2 = None
                if 'in/second_read' in run_ids_connections_files[run_id]:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]
                    input_fileset.append(r2)

                kallisto_eg = run.new_exec_group()

                kallisto = [self.get_tool('kallisto'), 'quant']
                kallisto.extend(['-i', str(self.get_option('i'))])
                out_path = run.get_output_directory_du_jour_placeholder()
                kallisto.extend(['-o', out_path])
                kallisto.extend(input_fileset)

                stderr_file = "%s-kallisto-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, input_fileset)
                stdout_file = "%s-kallisto-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, input_fileset)

                h5_file = "abundance.h5"
                run.add_output_file("abundance.h5", h5_file, input_fileset)
                tsv_file = "abundance.tsv"
                run.add_output_file("abundance.tsv", tsv_file, input_fileset)
                run_info = "run_info.json"
                run.add_output_file("run_info.json", run_info, input_fileset)

                kallisto_eg.add_command(kallisto, stdout_path=log_stdout,
                                        stderr_path=log_stderr)