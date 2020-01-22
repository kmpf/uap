from uaperrors import UAPError
import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PicardMarkDuplicates(AbstractStep):
    '''
    Identifies duplicate reads.
    This tool locates and tags duplicate reads (both PCR and optical/
    sequencing-driven) in a BAM or SAM file, where duplicate reads are defined
    as originating from the same original fragment of DNA.
    Duplicates are identified as read pairs having identical 5' positions
    (coordinate and strand) for both reads in a mate pair (and optinally,
    matching unique molecular identifier reads; see BARCODE_TAG option).
    Optical, or more broadly Sequencing, duplicates are duplicates that appear
    clustered together spatially during sequencing and can arise from optical/
    imagine-processing artifacts or from bio-chemical processes during clonal
    amplification and sequencing; they are identified using the READ_NAME_REGEX
    and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options.
    The tool's main output is a new SAM or BAM file in which duplicates have
    been identified in the SAM flags field, or optionally removed (see
    REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES), and optionally marked
    with a duplicate type in the 'DT' optional attribute.
    In addition, it also outputs a metrics file containing the numbers of
    READ_PAIRS_EXAMINED, UNMAPPED_READS, UNPAIRED_READS,
    UNPAIRED_READ_DUPLICATES, READ_PAIR_DUPLICATES, and
    READ_PAIR_OPTICAL_DUPLICATES.

    Usage example::

        java -jar picard.jar MarkDuplicates I=input.bam \
        O=marked_duplicates.bam M=marked_dup_metrics.txt

    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

    '''

    def __init__(self, pipeline):
        super(PicardMarkDuplicates, self).__init__(pipeline)

        self.set_cores(12)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/metrics')

        self.require_tool('picard-tools')

        # [Standard Picard Options:]

        self.add_option('TMP_DIR', str, optional = True,
                        description='A file. Default value: null. This option '
                        'may be specified 0 or more times.')
        self.add_option('VERBOSITY', str, optional=True,
                        choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                        description='Control verbosity of logging. '
                        'Default value: INFO. This option can be set to "null" '
                        'to clear the default value.')
        self.add_option('QUIET', bool, optional=True,
                        description='Whether to suppress job-summary info on '
                        'System.err. Default value: false. This option can be '
                        'set to "null" to clear the default value.')
        self.add_option('VALIDATION_STRINGENCY', str, optional=True,
                        choices=['STRICT', 'LENIENT', 'SILENT'],
                        description='Validation stringency for all SAM files '
                        'read by this program. Setting stringency to SILENT can '
                        'improve performance when processing a BAM file in '
                        'which variable-length data (read, qualities, tags) do '
                        'not otherwise need to be decoded. Default value: '
                        'STRICT. This option can be set to "null" to clear the '
                        'default value.')
        self.add_option('COMPRESSION_LEVEL', int, optional=True,
                        description='Compression level for all compressed files '
                        'created (e.g. BAM and GELI). Default value: 5. This '
                        'option can be set to "null" to clear the default value.')
        self.add_option('MAX_RECORDS_IN_RAM', int, optional=True,
                        description='When writing SAM files that need to be '
                        'sorted, this will specify the number of records stored '
                        'in RAM before spilling to disk. Increasing this number '
                        'reduces the number of file handles needed to sort a '
                        'SAM file, and increases the amount of RAM needed. '
                        'Default value: 500000. This option can be set to "null" '
                        'to clear the default value.')
        self.add_option('CREATE_INDEX', bool, optional=True,
                        description='Whether to create a BAM index when writing '
                        'a coordinate-sorted BAM file. Default value: false. '
                        'This option can be set to "null" to clear the default '
                        'value. ')
        self.add_option('CREATE_MD5_FILE', bool, optional=True,
                        description='Whether to create an MD5 digest for any '
                        'BAM or FASTQ files created. Default value: false. '
                        'This option can be set to "null" to clear the default '
                        'value.')
        self.add_option('REFERENCE_SEQUENCE', str, optional=True,
                        description='Reference sequence file. Default value: '
                        'null.')
        self.add_option('GA4GH_CLIENT_SECRETS', str, optional=True,
                        description='Google Genomics API client_secrets.json '
                        'file path. Default value: client_secrets.json. This '
                        'option can be set to "null" to clear the default '
                        'value.')

        # [Picard MarkDuplicates Options:]

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

        options = [
            # Standard Picard Options:
            'TMP_DIR', 'VERBOSITY', 'QUIET', 'VALIDATION_STRINGENCY',
            'COMPRESSION_LEVEL', 'MAX_RECORDS_IN_RAM', 'CREATE_INDEX',
            'CREATE_MD5_FILE', 'REFERENCE_SEQUENCE', 'GA4GH_CLIENT_SECRETS',
            # Picard MarkDuplicates Options:
            'PROGRAM_RECORD_ID', 'PROGRAM_GROUP_VERSION',
            'PROGRAM_GROUP_COMMAND_LINE', 'PROGRAM_GROUP_NAME',
            'COMMENT', 'ASSUME_SORTED', 'MAX_FILE_HANDLES',
            'SORTING_COLLECTION_SIZE_RATIO', 'READ_NAME_REGEX',
            'OPTICAL_DUPLICATE_PIXEL_DISTANCE'
        ]

        file_options = ['TMP_DIR', 'REFERENCE_SEQUENCE']

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
                value = str(self.get_option(option))
                if option in file_options:
                    value = os.path.abspath(value)
                option_list.append('%s=%s' % (option, value))

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/alignments']

                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise UAPError("Expected exactly one alignments file.")
                elif os.path.splitext(input_paths[0])[1] not in ['.sam', '.bam']:
                    raise UAPError(
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
                        mark_duplicates.extend(option_list)
                        exec_group.add_command(mark_duplicates)
