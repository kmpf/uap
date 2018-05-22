import sys
from logging import getLogger
from abstract_step import AbstractStep

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

        # output connections
        self.add_connection('out/fusion.txt')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        self.add_connection('out/tar_archive')

        # required tools
        self.require_tool('kallisto')
        self.require_tool('mkdir')
        self.require_tool('tar')
        self.require_tool('rm')

        # options
        self.add_option('cores', int, optional=True, default=6,
                        description="workaround to specify cores for grid \
                        engine and threads ie")

        self.add_option('index', str, optional=False, default=None,
                        description="Filename for the kallisto index to be \
                        used for quantification")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:

                # Get lst of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read'][0]
                sr_input = run_ids_connections_files[run_id]['in/second_read'][0]

                # Do we have paired end data and is it exactly one ?
                input_paths = [fr_input]

                if sr_input is None:
                    logger.error("Not paired end")
                    sys.exit(1)
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

                # kal_fusion = run.add_output_file(
                # 'fusion.txt', 'fusion.txt', input_paths)

                # Assemble kallisto command
                my_output = run.get_output_directory_du_jour_placeholder()
                with run.new_exec_group() as exec_group:

                    kallisto = [self.get_tool('kallisto'),
                                'quant',
                                '-i', self.get_option('index'),
                                '--fusion',
                                '-o', my_output,
                                fr_input,
                                sr_input]

                    exec_group.add_command(kallisto,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)
