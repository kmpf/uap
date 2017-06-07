from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class MergeGenecounts(AbstractStep):
    '''

    '''

    def __init__(self, pipeline):
        super(MergeGenecounts, self).__init__(pipeline)

        # in connections
        self.add_connection('in/first_read')

        # out connections
        self.add_connection('out/merged_counts')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                           engine and threads ie")

        self.add_option('t', int, optional=False,
                        description="tool name (htseq_count: htc, featureCounts: fc)")

        # required tools
        self.require_tool('merge_genecounts')

    def runs(self, run_ids_connections_files):

        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                print(run)