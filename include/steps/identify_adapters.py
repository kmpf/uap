from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class IdentifyAdapters(AbstractStep):
    '''
    Uses AdapterRemoval to identify adapter sequences from paired read data.

    AdapterRemoval (ver. 2.1.7)
    This program searches for and removes remnant adapter sequences from
    your read data.  The program can analyze both single end and paired end
    data.  For detailed explanation of the parameters, please refer to the
    man page.  For comments, suggestions  and feedback please contact Stinus
    Lindgreen (stinus@binf.ku.dk) and Mikkel Schubert (MikkelSch@gmail.com).
    If you use the program, please cite the paper:
    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid
    adapter trimming, identification, and read merging.
    BMC Research Notes, 12;9(1):88.
    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2
    "Pipeline specific "input and output expected to be gzipped"
    '''

    def __init__(self, pipeline):
        super(IdentifyAdapters, self).__init__(pipeline)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.require_tool('adapterremoval')

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                r2 = run_ids_connections_files[run_id]['in/second_read'][0]

                ar_exec_group = run.new_exec_group()
                ar = [self.get_tool('adapterremoval')]
                ar.append('--identify-adapters')
                ar.extend(['--file1', r1, '--file2', r2])

                stderr_file = "%s-adapterremoval-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, [r1, r2])
                stdout_file = "%s-adapterremoval-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, [r1, r2])

                ar_exec_group.add_command(ar, stdout_path=log_stdout,
                                          stderr_path=log_stderr)
