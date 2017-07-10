from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Salmon(AbstractStep):
    '''
    Salmon


    '''

    def __init__(self, pipeline):
        super(Salmon, self).__init__(pipeline)

        # input connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        # output connections
        self.add_connection('out/cmd_info.json')
        self.add_connection('out/lib_format_counts.json')
        self.add_connection('out/quant.sf')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        self.add_connection('out/flenDist.txt')
        #self.add_connection('out/salmon_quant.log')
        #self.add_connection('out/logs')

        # required tools
        self.require_tool('salmon')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                        engine and threads ie")

        self.add_option('i', str, optional=False, default=None,
                        description="Salmon index")

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

                salmon_eg = run.new_exec_group()

                salmon = [self.get_tool('salmon'), 'quant']
                salmon.extend(['-i', str(self.get_option('i'))])
                salmon.extend(['-l', 'ISF'])
                out_path = run.get_output_directory_du_jour_placeholder()
                salmon.extend(['-o', out_path])
                salmon.extend(['-1', r1])
                (r2 is not None) and (salmon.extend(['-2', r2]))

                stderr_file = "%s-salmon-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, input_fileset)
                stdout_file = "%s-salmon-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, input_fileset)

                run.add_output_file("cmd_info.json", "cmd_info.json", input_fileset)
                run.add_output_file("lib_format_counts.json", "lib_format_counts.json", input_fileset)
                run.add_output_file("quant.sf", "quant.sf", input_fileset)

                run.add_output_file("flenDist.txt", "libParams/flenDist.txt", input_fileset)
                #run.add_output_file("salmon_quant.log", "salmon_quant.log", input_fileset)
                #run.add_output_file("logs", "logs", input_fileset)

                salmon_eg.add_command(salmon, stdout_path=log_stdout,
                                        stderr_path=log_stderr)