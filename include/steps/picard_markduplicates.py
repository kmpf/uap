from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')

class PicardMarkDuplicates(AbstractStep):
    '''
    Identifies duplicate reads.

    This tool locates and tags duplicate reads in a BAM or SAM file, where
    duplicate reads are defined as originating from a single fragment of DNA.
    Duplicates can arise during sample preparation e.g. library construction
    using PCR. See also EstimateLibraryComplexity for additional notes on PCR
    duplication artifacts. Duplicate reads can also result from a single
    amplification cluster, incorrectly detected as multiple clusters by the
    optical sensor of the sequencing instrument. These duplication artifacts
    are referred to as optical duplicates.

    The MarkDuplicates tool works by comparing sequences in the 5 prime
    positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG
    option is available to facilitate duplicate marking using molecular
    barcodes. After duplicate reads are collected, the tool differentiates the
    primary and duplicate reads using an algorithm that ranks reads by the sums
    of their base-quality scores (default method).

    The tool's main output is a new SAM or BAM file, in which duplicates have
    been identified in the SAM flags field for each read. Duplicates are marked
    with the hexadecimal value of 0x0400, which corresponds to a decimal value
    of 1024. If you are not familiar with this type of annotation, please see
    the following blog post for additional information.

    Although the bitwise flag annotation indicates whether a read was marked as
    a duplicate, it does not identify the type of duplicate. To do this, a new
    tag called the duplicate type (DT) tag was recently added as an optional
    output in the 'optional field' section of a SAM/BAM file. Invoking the
    TAGGING_POLICY option, you can instruct the program to mark all the
    duplicates (All), only the optical duplicates (OpticalOnly), or no
    duplicates (DontTag). The records within the output of a SAM/BAM file will
    have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as
    either library/PCR-generated duplicates (LB), or sequencing-platform
    artifact duplicates (SQ). This tool uses the READ_NAME_REGEX and the
    OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify
    and differentiate duplicate types. Set READ_NAME_REGEX to null to skip
    optical duplicate detection, e.g. for RNA-seq or other data where duplicate
    sets are extremely large and estimating library complexity is not an aim.
    Note that without optical duplicate counts, library size estimation will be
    inaccurate.

    MarkDuplicates also produces a metrics file indicating the numbers of
    duplicates for both single- and paired-end reads.

    The program can take either coordinate-sorted or query-sorted inputs,
    however the behavior is slightly different. When the input is
    coordinate-sorted, unmapped mates of mapped records and
    supplementary/secondary alignments are not marked as duplicates. However,
    when the input is query-sorted (actually query-grouped), then unmapped
    mates and secondary/supplementary reads are not excluded from the
    duplication test and can be marked as duplicate reads.

    If desired, duplicates can be removed using the REMOVE_DUPLICATE and
    REMOVE_SEQUENCING_DUPLICATES options.

    Usage example::

        java -jar picard.jar MarkDuplicates \
             I=input.bam \
             O=marked_duplicates.bam \
             M=marked_dup_metrics.txt

    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

    '''

    def __init__(self, pipeline):
        super(PicardMarkDuplicates, self).__init__(pipeline)

        self.set_cores(12)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/metrics')

        # Step was tested for picard-tools release 1.113
        self.require_tool('picard-tools')

        # [Standard Picard Options:]

        self.add_option('TMP_DIR', str, optional = True,
                        description = 'A file. Default value: null. This option '
                        'may be specified 0 or more times.')
        self.add_option('VERBOSITY', str, optional = True,
                        choices = ['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                        description = 'Control verbosity of logging. '
                        'Default value: INFO. This option can be set to "null" '
                        'to clear the default value.')
        self.add_option('QUIET', bool, optional = True,
                        description = 'Whether to suppress job-summary info on '
                        'System.err. Default value: false. This option can be '
                        'set to "null" to clear the default value.')
        self.add_option('VALIDATION_STRINGENCY', str, optional=True,
                        choices=['STRICT', 'LENIENT', 'SILENT', 'null'],
                        description='Validation stringency for all SAM files '
                        'read by this program. Setting stringency to SILENT can '
                        'improve performance when processing a BAM file in '
                        'which variable-length data (read, qualities, tags) do '
                        'not otherwise need to be decoded. Default value: '
                        'STRICT. This option can be set to "null" to clear the '
                        'default value.')
        self.add_option('COMPRESSION_LEVEL', int, optional=True,
                        description='Compression level for all compressed files '
                        'created (e.g. BAM and GELI). Default value: 5.')
        self.add_option('MAX_RECORDS_IN_RAM', int, optional=True,
                        description='When writing SAM files that need to be '
                        'sorted, this will specify the number of records stored '
                        'in RAM before spilling to disk. Increasing this number '
                        'reduces the number of file handles needed to sort a '
                        'SAM file, and increases the amount of RAM needed. '
                        'Default value: 500000.')
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

        self.add_option('BARCODE_TAG', str, optional = True,
                        description = "Barcode SAM tag (ex. BC for 10X "
                        "Genomics) Default value: null.")
        self.add_option('READ_ONE_BARCODE_TAG', str, optional = True,
                        description = "Read one barcode SAM tag (ex. BX for "
                        "10X Genomics) Default value: null.")
        self.add_option('READ_TWO_BARCODE_TAG', str, optional = True,
                        description = "Read two barcode SAM tag (ex. BX for "
                        "10X Genomics) Default value: null.")
        self.add_option('TAG_DUPLICATE_SET_MEMBERS', str, optional = True,
                        choices = ["true", "false", "null"], description = "If "
                        "a read appears in a duplicate set, add two tags. The "
                        "first tag, DUPLICATE_SET_SIZE_TAG (DS), indicates the "
                        "size of the duplicate set. The smallest possible DS "
                        "value is 2 which occurs when two reads map to the "
                        "same portion of the reference only one of which is "
                        "marked as duplicate. The second tag, "
                        "DUPLICATE_SET_INDEX_TAG (DI), represents a unique "
                        "identifier for the duplicate set to which the record "
                        "belongs. This identifier is the index-in-file of the "
                        "representative read that was selected out of the "
                        "duplicate set. Default value: false. This option can "
                        "be set to 'null' to clear the default value. Possible "
                        "values: {true, false}")
        self.add_option('REMOVE_SEQUENCING_DUPLICATES', str, optional = True,
                        choices = ["true", "false", "null"], description = "If "
                        "true remove 'optical' duplicates and other duplicates "
                        "that appear to have arisen from the sequencing "
                        "process instead of the library preparation process, "
                        "even if REMOVE_DUPLICATES is false. If "
                        "REMOVE_DUPLICATES is true, all duplicates are removed "
                        "and this option is ignored. Default value: false. "
                        "This option can be set to 'null' to clear the default "
                        "value. Possible values: {true, false}")
        self.add_option('TAGGING_POLICY', str, optional = True,
                        choices = ["DontTag", "OpticalOnly", "All"],
                        description = "Determines how duplicate types are "
                        "recorded in the DT optional attribute. Default value: "
                        "DontTag. This option can be set to 'null' to clear "
                        "the default value. Possible values: {DontTag, "
                        "OpticalOnly, All}")
        self.add_option('ASSUME_SORT_ORDER', str, optional = True,
                        choices = ["unsorted", "queryname", "coordinate",
                                   "duplicate", "unknown"],
                        description = "If not null, assume that the input file "
                        "has this order even if the header says otherwise. "
                        "Default value: null. Possible values: {unsorted, "
                        "queryname, coordinate, duplicate, unknown}")
        self.add_option('DUPLICATE_SCORING_STRATEGY', str, optional = True,
                        choices = ["SUM_OF_BASE_QUALITIES",
                                   "TOTAL_MAPPED_REFERENCE_LENGTH", "RANDOM"],
                        description = "The scoring strategy for choosing the "
                        "non-duplicate among candidates. This option can be "
                        "set to 'null' to clear the default value. Default "
                        "value: SUM_OF_BASE_QUALITIES. Possible values: "
                        "{SUM_OF_BASE_QUALITIES, "
                        "TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}")

        # Options supported by 'picard-tools MarkDuplicates' release 1.113
        self.add_option('PROGRAM_RECORD_ID', str, optional = True, description =
                        "The program record ID for the @PG record(s) created by "
                        "this program. Set to null to disable PG record "
                        "creation. This string may have a suffix appended to "
                        "avoid collision with other program record IDs. "
                        "Default value: MarkDuplicates. This option can be set "
                        "to 'null' to clear the default value.")
        self.add_option('PROGRAM_GROUP_VERSION', str, optional = True,
                        description = "Value of VN tag of PG record to be "
                        "created. If not specified, the version will be "
                        "detected automatically. Default value: null.")
        self.add_option('PROGRAM_GROUP_COMMAND_LINE', str, optional = True,
                        description = "Value of CL tag of PG record to be "
                        "created. If not supplied the command line will be "
                        "detected automatically. Default value: null.")
        self.add_option('PROGRAM_GROUP_NAME', str, optional = True,
                        description = "Value of PN tag of PG record to be "
                        "created. Default value: MarkDuplicates. This option "
                        "can be set to 'null' to clear the default value.")
        self.add_option('COMMENT', str, optional = True,
                        description = "Comment(s) to include in the output "
                        "file's header. Default value: null. This option may "
                        "be specified 0 or more times.")
        self.add_option('MAX_FILE_HANDLES', int, optional = True)
        self.add_option('REMOVE_DUPLICATES', str, optional = True,
                        default = "true",
                        choices = ["true", "null", "false"],
                        description = "If true do not write duplicates to the "
                        "output file instead of writing them with appropriate "
                        "flags set. Default value: false. This option can be "
                        "set to 'null' to clear the default value. Possible "
                        "values: {true, false}")
        self.add_option('ASSUME_SORTED', str, optional = True,
                        choices = ["true", "null", "false"],
                        description = "If true, assume that the input file is "
                        "coordinate sorted even if the header says otherwise. "
                        "Default value: false. This option can be set to 'null' "
                        "to clear the default value. Possible values: "
                        "{true, false}")
        self.add_option('MAX_FILE_HANDLES_FOR_READ_ENDS_MAP', int,
                        optional = True,
                        description = "Maximum number of file "
                        "handles to keep open when spilling read ends to disk. "
                        "Set this number a little lower than the per-process "
                        "maximum number of file that may be open. This number "
                        "can be found by executing the 'ulimit -n' command on "
                        "a Unix system. Default value: 8000.")
        self.add_option('SORTING_COLLECTION_SIZE_RATIO', float,
                        optional = True,
                        description = "This number, plus the maximum RAM "
                        "available to the JVM, determine the memory footprint "
                        "used by some of the sorting collections. If you are "
                        "running out of memory, try reducing this number. "
                        "Default value: 0.25.")
        self.add_option('READ_NAME_REGEX', str, optional = True,
                        description = "Regular expression that can be used to "
                        "parse read names in the incoming SAM file. Read names "
                        "are parsed to extract three variables: tile/region, x "
                        "coordinate and y coordinate. These values are used to "
                        "estimate the rate of optical duplication in order to "
                        "give a more accurate estimated library size. Set this "
                        "option to null to disable optical duplicate detection, "
                        "e.g. for RNA-seq or other data where duplicate sets "
                        "are extremely large and estimating library complexity "
                        "is not an aim. Note that without optical duplicate "
                        "counts, library size estimation will be inaccurate. "
                        "The regular expression should contain three capture "
                        "groups for the three variables, in order. It must "
                        "match the entire read name. Note that if the default "
                        "regex is specified, a regex match is not actually "
                        "done, but instead the read name is split on colon "
                        "character. For 5 element names, the 3rd, 4th and 5th "
                        "elements are assumed to be tile, x and y values. For "
                        "7 element names (CASAVA 1.8), the 5th, 6th, and 7th "
                        "elements are assumed to be tile, x and y values. "
                        "Default value: . This option can be set to 'null' to "
                        "clear the default value.")
        self.add_option('OPTICAL_DUPLICATE_PIXEL_DISTANCE', int, optional = True,
                        description = "The maximum offset between two duplicate "
                        "clusters in order to consider them optical duplicates. "
                        "The default is appropriate for unpatterned versions of "
                        "the Illumina platform. For the patterned flowcell "
                        "models, 2500 is more appropriate. For other platforms "
                        "and models, users should experiment to find what works "
                        "best. Default value: 100. This option can be set to "
                        "'null' to clear the default value.")

    def runs(self, run_ids_connections_files):

        options = [
            # Standard Picard Options:
            'TMP_DIR', 'VERBOSITY', 'QUIET', 'VALIDATION_STRINGENCY',
            'COMPRESSION_LEVEL', 'MAX_RECORDS_IN_RAM', 'CREATE_INDEX',
            'CREATE_MD5_FILE', 'REFERENCE_SEQUENCE', 'GA4GH_CLIENT_SECRETS',
            # Picard MarkDuplicates Options:
            'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP',
            'SORTING_COLLECTION_SIZE_RATIO', 'BARCODE_TAG',
            'READ_ONE_BARCODE_TAG', 'READ_TWO_BARCODE_TAG',
            'TAG_DUPLICATE_SET_MEMBERS', 'REMOVE_SEQUENCING_DUPLICATES',
            'TAGGING_POLICY', 'REMOVE_DUPLICATES', 'ASSUME_SORT_ORDER',
            'DUPLICATE_SCORING_STRATEGY',
            'PROGRAM_RECORD_ID', 'PROGRAM_GROUP_VERSION',
            'PROGRAM_GROUP_COMMAND_LINE', 'PROGRAM_GROUP_NAME',
            'COMMENT', 'ASSUME_SORTED',
            'READ_NAME_REGEX', 'OPTICAL_DUPLICATE_PIXEL_DISTANCE'
        ]

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
                    raise StepError(self, "Expected exactly one alignments file.")
                elif os.path.splitext(input_paths[0])[1] not in ['.sam', '.bam']:
                    raise StepError(self,
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
                            'METRICS_FILE=%s' % metrics
                        ]
                        mark_duplicates.extend(option_list)
                        exec_group.add_command(mark_duplicates)
