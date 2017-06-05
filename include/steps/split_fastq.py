from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class SplitFastq(AbstractStep):
    '''

    '''

    def __init__(self, pipeline):
        super(SplitFastq, self).__init__(pipeline)

        # in connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        # out connections
        self.add_connection('out/out_first_read')
        # self.add_connection('out/out_second_read')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                    engine and threads ie")

        self.add_option('readcount', int, optional=False,
                        description="Number of reads per targetfile")

        # required tools
        self.require_tool('split_fastq')

    def runs(self, run_ids_connections_files):

        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():

            # TODO: create new run id

            with self.declare_run(run_id) as run:
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                readcount = self.get_option('readcount')
                print(r1, readcount)
                output_fileset = [r1]

                split_fastq = [self.get_tool('split_fastq'), '-i', r1, '-n', readcount]

                if run_ids_connections_files[run_id]['in/second_read'][0] == None:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]