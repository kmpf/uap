from uaperrors import UAPError
from logging import getLogger
from abstract_step import AbstractStep
import os

logger = getLogger('uap_logger')


class Tcount2gcount(AbstractStep):
    '''

    '''

    def __init__(self, pipeline):
        super(Tcount2gcount, self).__init__(pipeline)

        # in connections
        self.add_connection('in/counts')
        self.add_connection('in/annotation')

        # out connections
        self.add_connection('out/counts')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                           engine and threads ie")

        self.add_option('t', str, optional=True,
                        description='source tool of input file. ' \
                                    'Possible values: kallisto, salmon')

        self.add_option('m', str, optional=True,
                        description='transcript to gene mapping file. ' \
                        'Required Format example (per row): ' \
                        'ENST00000527779.1\tENSG00000137692.11' \
                        ' or gtf')


        self.add_option('kallisto-extended', bool, optional=True,
                        description='writes extended format includign tpm. ' )


        # required tools
        self.require_tool('tcount2gcount')

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))


        annotation = None
        for run_id in run_ids_connections_files.keys():
            if 'in/annotation' in run_ids_connections_files[run_id]:
                annotation = run_ids_connections_files[run_id]['in/annotation'][0]

        for run_id in run_ids_connections_files.keys():
            if 'in/annotation' in run_ids_connections_files[run_id]:
                continue 

            with self.declare_run(run_id) as run:
                counts = run_ids_connections_files[run_id]['in/counts'][0]

                tool_name = self.get_option('t')
                file_name = '%s-gene-abundance.tsv' % (run_id)
                run.add_output_file('counts', file_name, [counts])
                basename = run.get_output_directory_du_jour_placeholder()

                cmd = [self.get_tool('tcount2gcount')]

                if self.is_option_set_in_config('m'): 
                     cmd.extend(['-m', os.path.abspath(self.get_option('m'))])
                else:                    
                    if annotation:
                        cmd.extend(['-m', os.path.abspath(annotation)])
                    else:  
                        raise UAPError("%s no annotation give via config or connection" % run_id) 
                    
                if self.is_option_set_in_config('kallisto-extended'): 
                    cmd.append('--kallisto-extended')

                cmd.extend(['-i', counts,
                            '-t', tool_name,
                            '-o', basename + '/' + file_name])

                convert_exec_group = run.new_exec_group()
                convert_exec_group.add_command(cmd)
