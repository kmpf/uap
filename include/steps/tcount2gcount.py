from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Tcount2gcount(AbstractStep):
    '''

    '''

    def __init__(self, pipeline):
        super(Tcount2gcount, self).__init__(pipeline)

        # in connections
        self.add_connection('in/transcript_counts')

        # out connections
        self.add_connection('out/gene_counts')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                           engine and threads ie")

        self.add_option('t', str, optional=False,
                        description='source tool of input file. ' \
                                    'Possible values: kallisto, salmon')

        self.add_option('m', str, optional=False,
                        description='transcript to gene mapping file. ' \
                        'Required Format example (per row): ' \
                        'ENST00000527779.1\tENSG00000137692.11')

        # required tools
        self.require_tool('tcount2gcount')

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        t_counts = []
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                t_counts.append(run_ids_connections_files[run_id]['in/transcript_counts'][0])

        for i in range(0, len(t_counts)):
            new_run_id = 'tcount2gcount_%s' % i
            with self.declare_run(new_run_id) as run:
                tool_name = self.get_option('t')
                file_name = '%s_%s.tsv' % (new_run_id, tool_name)
                run.add_output_file('gene_counts', file_name, t_counts)
                basename = run.get_output_directory_du_jour_placeholder()
                convert_command = [
                    self.get_tool('tcount2gcount'),
                    '-m', self.get_option('m'),
                    '-i', t_counts[i],
                    '-t', tool_name,
                    '-o', basename + '/' + file_name
                ]
                convert_exec_group = run.new_exec_group()
                convert_exec_group.add_command(convert_command)