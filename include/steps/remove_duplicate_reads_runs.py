import sys
from abstract_step import *
import misc
import process_pool
import yaml


class RemoveDuplicates(AbstractStep):
    '''
    Duplicates are removed by Picard tools 'MarkDuplicates'.

    typical command line::

        MarkDuplicates INPUT=<SAM/BAM> OUTPUT=<SAM/BAM>
                       METRICS_FILE=<metrics-out> REMOVE_DUPLICATES=true
    '''

    def __init__(self, pipeline):
        super(RemoveDuplicates, self).__init__(pipeline)

        self.set_cores(12)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/metrics')

        self.require_tool('MarkDuplicates')

    def runs(self, run_ids_connections_files):
        # Compile list of options
        options = []
        # Compile list of set options
        set_options = [option for option in options if
                       self.is_option_set_in_config(option)]

        for run_id in run_ids_connections_files.keys():
            # Get input alignments
            input_paths = run_ids_connections_files[run_id]['in/alignments']
            with self.declare_run(run_id) as run:
                for input_path in input_paths:
                    with run.new_exec_group() as exec_group:
                        # Assemble MarkDuplicates command
                        mark_duplicates = [
                            self.get_tool('MarkDuplicates'),
                            'INPUT=%s' % input_path,
                            'OUTPUT=%s' % run.add_output_file(
                                'alignments',
                                '%s-rm-dup.bam' % run_id,
                                input_path),
                            'METRICS_FILE=%s' % run.add_output_file(
                                'metrics',
                                '%s-rm-dup-metrics.txt' % run_id,
                                input_path),
                            'REMOVE_DUPLICATES=true'
                        ]
                        exec_group.add_command(
                            mark_duplicates,
                        )
