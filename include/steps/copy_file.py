import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class CopyFile(AbstractStep):
    '''
    copies a file or a list of files defined by there
    dependencies and filenames
    '''

    def __init__(self, pipeline):
        super(CopyFile, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/sequence')
        self.add_connection('out/copied')

        self.require_tool('cp')

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/sequence']

                with run.new_exec_group() as cp_exec_group:
                    for input_file in input_paths:
                        file_name = os.path.basename(input_file)
                        out_file = run.add_output_file('copied',
                                                       file_name,
                                                       input_paths)
                        cp = [self.get_tool('cp'),
                              input_file, out_file]
                        cp_exec_group.add_command(cp)
