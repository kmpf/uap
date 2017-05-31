from Bio import SeqIO

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
        self.add_connection('out/out_second_read')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                    engine and threads ie")

        self.add_option('readcount', int, optional=False,
                        description="Number of reads per targetfile")

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                print('Test')