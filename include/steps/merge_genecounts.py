from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class MergeGenecounts(AbstractStep):
    '''

    '''

    def __init__(self, pipeline):
        super(MergeGenecounts, self).__init__(pipeline)

        # in connections
        self.add_connection('in/counts')

        # out connections
        self.add_connection('out/merged_counts')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                           engine and threads ie")

        self.add_option('t', str, optional=False,
                        description="tool name (htseq_count: htc, featureCounts: fc)")

        # required tools
        self.require_tool('merge_genecounts')

    def runs(self, run_ids_connections_files):

        self.set_cores(self.get_option('cores'))

        gc_files = []
        for run_id in run_ids_connections_files.keys():
            if run_ids_connections_files[run_id]['in/counts']:
                with self.declare_run(run_id) as run:
                    # collect needed files
                    gc_file = run_ids_connections_files[run_id]['in/counts'][0]
                    gc_files.append(gc_file)

        new_run_id = 'merged_genecounts'
        tool_name = self.get_option('t')
        with self.declare_run(new_run_id) as run:
            file_name = '%s_%s.txt' % (new_run_id, tool_name)
            run.add_output_file('merged_counts', file_name, gc_files)
            merge_command = [
                self.get_tool('merge_genecounts'),
                '-t', tool_name,
                '-o', '.',
                '-p', file_name,
                ' '.join(gc_files)
            ]
            merge_exec_group = run.new_exec_group()
            merge_exec_group.add_command(merge_command)
