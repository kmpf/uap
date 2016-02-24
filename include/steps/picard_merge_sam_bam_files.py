import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PicardMergeSamFiles(AbstractStep):
    '''

    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles

    '''

    def __init__(self, pipeline):
        super(PicardMergeSamFiles, self).__init__(pipeline)

        self.set_cores(12)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('ln')
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

        # [Picard MergeSamFiles Options:]

        self.add_option('SORT_ORDER', str, optional = True,
                        choices=['unsorted', 'queryname', 'coordinate',
                                 'duplicate'],
                        description="Sort order of output file. Default value: "
                        "coordinate. This option can be set to 'null' to clear "
                        "the default value. Possible values: {unsorted, "
                        "queryname, coordinate, duplicate}")
        self.add_option('ASSUME_SORTED', bool, optional=True,
                        description="If true, assume that the input files are "
                        "in the same sort order as the requested output sort "
                        "order, even if their headers say otherwise. Default "
                        "value: false. This option can be set to 'null' to "
                        "clear the default value. Possible values: "
                        "{true, false}")
        self.add_option('MERGE_SEQUENCE_DICTIONARIES', bool, optional=True,
                        description="Merge the sequence dictionaries. Default "
                        "value: false. This option can be set to 'null' to "
                        "clear the default value. Possible values: "
                        "{true, false}")
        self.add_option('USE_THREADING', bool, optional=True,
                        description="Option to create a background thread to "
                        "encode, compress and write to disk the output file. "
                        "The threaded version uses about 20% more CPU and "
                        "decreases runtime by ~20% when writing out a "
                        "compressed BAM file. Default value: false. This "
                        "option can be set to 'null' to clear the default "
                        "value. Possible values: {true, false}")
        self.add_option('COMMENT', str, optional=True,
                        description="Comment(s) to include in the merged "
                        "output file's header. Default value: null.")
        self.add_option('INTERVALS', str, optional=True,
                        description="An interval list file that contains the "
                        "locations of the positions to merge. Assume bam are "
                        "sorted and indexed. The resulting file will contain "
                        "alignments that may overlap with genomic regions "
                        "outside the requested region. Unmapped reads are "
                        "discarded. Default value: null.")

    def runs(self, run_ids_connections_files):

        options = [
            # Standard Picard Options:
            'TMP_DIR', 'VERBOSITY', 'QUIET', 'VALIDATION_STRINGENCY',
            'COMPRESSION_LEVEL', 'MAX_RECORDS_IN_RAM', 'CREATE_INDEX',
            'CREATE_MD5_FILE', 'REFERENCE_SEQUENCE', 'GA4GH_CLIENT_SECRETS',
            # Picard MarkDuplicates Options:
            'SORT_ORDER', 'ASSUME_SORTED', 'MERGE_SEQUENCE_DICTIONARIES',
            'USE_THREADING', 'COMMENT', 'INTERVALS'
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
                elif os.path.splitext(input_paths[0])[1] not in ['.sam', '.bam']:
                    logger.error(
                        "The file %s seems not to be a SAM or BAM file. At "
                        "least the suffix is wrong." % input_paths[0]
                    )
                    sys.exit(1)
                elif self.is_option_set_in_config("INTERVALS") and \
                     not os.path.exists(self.get_option("INTERVALS")):
                    logger.error("The path %s given to option 'INTERVALS' is "
                                 "not pointing to a file.")
                    sys.exit(1)
                elif len(input_paths) == 0:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) == 1:
                    base = os.path.basename(input_paths[0])
                    with run.new_exec_group() as ln_alignment:
                        # 1. command: Create symbolic link to original bam file
                        # (use absolute path)
                        ln = [self.get_tool('ln'), '-s',
                              input_paths[0],
                              run.add_output_file(
                                  'alignments', base, input_paths)]
                        ln_alignment.add_command(ln)

                else:
                    with run.new_exec_group() as exec_group:
                        alignments = run.add_output_file(
                            'alignments', '%s-merged.bam' % run_id,
                            input_paths)
                        merge_sam_files = [
                            self.get_tool('picard-tools'),
                            'MergeSamFiles']
                        for f in input_paths:
                            merge_sam_files.append('INPUT=%s' % f)
                        merge_sam_files.append('OUTPUT=%s' % alignments)
                        merge_sam_files.extend(option_list)
                        exec_group.add_command(merge_sam_files)
