import sys
from abstract_step import *
import misc
import process_pool
import yaml


class PicardMarkDuplicates(AbstractStep):
    '''
    Documentation:
    http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates

    Examines aligned records in the supplied SAM or BAM file to locate duplicate
    molecules. All records are then written to the output file with the
    duplicate records flagged.
    '''

    def __init__(self, pipeline):
        super(PicardMarkDuplicates, self).__init__(pipeline)

        self.set_cores(12)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/metrics')
        
        self.require_tool('picard-tools')

        self.add_option('PROGRAM_RECORD_ID', str, optional = True)
        self.add_option('PROGRAM_GROUP_VERSION', str, optional = True)
        self.add_option('PROGRAM_GROUP_COMMAND_LINE', str, optional = True)
        self.add_option('PROGRAM_GROUP_NAME', str, optional = True)
        self.add_option('COMMENT', str, optional = True)

        self.add_option('ASSUME_SORTED', bool, optional = True)
        self.add_option('MAX_FILE_HANDLES', int, optional = True)
        self.add_option('SORTING_COLLECTION_SIZE_RATIO', float, optional = True)
        self.add_option('READ_NAME_REGEX', str, optional = True)
        self.add_option('OPTICAL_DUPLICATE_PIXEL_DISTANCE', int, optional = True)

    def runs(self, run_ids_connections_files):

        options = ['PROGRAM_RECORD_ID', 'PROGRAM_GROUP_VERSION',
                   'PROGRAM_GROUP_COMMAND_LINE', 'PROGRAM_GROUP_NAME',
                   'COMMENT', 'ASSUME_SORTED', 'MAX_FILE_HANDLES',
                   'SORTING_COLLECTION_SIZE_RATIO', 'READ_NAME_REGEX',
                   'OPTICAL_DUPLICATE_PIXEL_DISTANCE']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('%s=true' % option)
                else:
                    option_list.append('%s=false' % option)
            else:
                option_list.append(
                    '%s=%s' % (option, str(self.get_option(option)))
                )

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/alignments']

                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                elif os.path.splitext(input_paths[0])[1] not in ['.sam', '.bam']:
                    raise StandardError(
                        "The file %s seems not to be a SAM or BAM file. At "
                        "least the suffix is wrong." % input_paths[0]
                    )
                else:
                    with run.new_exec_group() as exec_group:
                        alignments = run.add_output_file(
                            'alignments', '%s-rm-dup.bam' % run_id, input_paths)
                        metrics = run.add_output_file(
                            "metrics", '%s-rm-dup-metrics.txt' % run_id,
                            input_paths)
                        mark_duplicates = [
                            self.get_tool('picard-tools'), 'MarkDuplicates',
                            'INPUT=%s' % input_paths[0],
                            'OUTPUT=%s' % alignments,
                            'METRICS_FILE=%s' % metrics,
                            'REMOVE_DUPLICATES=true'                
                        ]
                        exec_group.add_command(mark_duplicates)
