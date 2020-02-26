import logging
from abstract_step import AbstractStep
import os
import sys

logger = logging.getLogger('uap_logger')


class CollectScs(AbstractStep):
    '''
    supa custom succesive aligment info gather thingy
    '''

    def __init__(self, pipeline):
        super(CollectScs, self).__init__(pipeline)

        self.set_cores(1)  # muss auch in den Decorator
        self.add_connection('in/scs_metrics')
        self.add_connection('out/yaml')

        self.require_tool('collect_scs')

        self.add_option('types', list, default=None, help ="genome names of alignments")

        self.add_option('rrna-aln-pos', str,   default=None,
                        help ="rRNA counts position in succesive alignment", optional=False)

        self.add_option('library-type', str,
                        choices=['unstranded', 'firststranded', 'secondstranded'],
                        default=None, help ="See tophat manual.", optional=False)


    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/scs_metrics']
                if input_paths == [None]:
                    run.add_empty_output_connection("yaml")

                out = run.add_output_file(
                    "yaml",
                    "%s-scs.yaml" %  (run_id),
                    input_paths)

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
                        scs = [self.get_tool('collect_scs'),
                               '--name', run_id,
                               '--outfile', out,
                               '--rrna-aln-pos', self.get_option('rrna-aln-pos'),
                               '--library-type', self.get_option('library-type'),
                               '--types']
                        scs.extend(self.get_option('types'))
                        scs.append('--infiles')
                        scs.extend(input_paths)

                        pipe.add_command(scs)


