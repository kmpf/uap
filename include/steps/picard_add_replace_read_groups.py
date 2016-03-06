import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PicardAddOrReplaceReadGroups(AbstractStep):
    '''
    Replace read groups in a BAM file. This tool enables the user to replace all
    read groups in the INPUT file with a single new read group and assign all
    reads to this read group in the OUTPUT BAM file.

    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups

    '''

    def __init__(self, pipeline):
        super(PicardAddOrReplaceReadGroups, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
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

        # [Picard AddOrReplaceReadGroups Options:]

        self.add_option('SORT_ORDER', str, optional = True,
                        choices=['unsorted', 'queryname', 'coordinate',
                                 'duplicate'],
                        description="Optional sort order to output in. If not "
                        "supplied OUTPUT is in the same order as INPUT. "
                        "Default value: null. Possible values: {unsorted, "
                        "queryname, coordinate, duplicate}")
        self.add_option('RGID', str, optional=True, description="Read Group "
                        "ID Default value: 1. This option can be set to 'null' "
                        "to clear the default value.")
        self.add_option('RGLB', str, optional=False,
                        description="Read Group library")
        self.add_option('RGPL', str, optional=False,
                        description="Read Group platform (e.g. illumina, solid)"
        )
        self.add_option('RGPU', str, optional=False,
                        description="Read Group platform unit (eg. run barcode)"
        )
        self.add_option('RGCN', str, optional=True, description="Read Group "
                        "sequencing center name. Default value: null.")
        self.add_option('RGDS', str, optional=True, description="Read Group "
                        "description. Default value: null.")
        self.add_option('RGDT', str, optional=True, description="Read Group "
                        "run date. Default value: null.")
        self.add_option('RGPI', int, optional=True, description="Read Group "
                        "predicted insert size. Default value: null.")
        self.add_option('RGPG', str, optional=True, description="Read Group "
                        "program group. Default value: null.")
        self.add_option('RGPM', str, optional=True, description="Read Group "
                        "platform model. Default value: null.")

    def runs(self, run_ids_connections_files):

        options = [
            # Standard Picard Options:
            'TMP_DIR', 'VERBOSITY', 'QUIET', 'VALIDATION_STRINGENCY',
            'COMPRESSION_LEVEL', 'MAX_RECORDS_IN_RAM', 'CREATE_INDEX',
            'CREATE_MD5_FILE', 'REFERENCE_SEQUENCE', 'GA4GH_CLIENT_SECRETS',
            # Picard MarkDuplicates Options:
            'SORT_ORDER', 'RGID', 'RGLB', 'RGPL', 'RGPU', 'RGCN',
            'RGDS', 'RGDT', 'RGPI', 'RGPG', 'RGPM'
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
                    logger.error("Expected exactly one alignments file.")
                    sys.exit(1)
                elif os.path.splitext(input_paths[0])[1] not in ['.sam', '.bam']:
                    logger.error(
                        "The file %s seems not to be a SAM or BAM file. At "
                        "least the suffix is wrong." % input_paths[0]
                    )
                    sys.exit(1)
                else:
                    with run.new_exec_group() as exec_group:
                        alignments = run.add_output_file(
                            'alignments', os.path.basename(input_paths[0]),
                            input_paths)
                        add_replace_read_groups = [
                            self.get_tool('picard-tools'),
                            'AddOrReplaceReadGroups',
                            'INPUT=%s' % input_paths[0],
                            'OUTPUT=%s' % alignments,
                            'RGSM=%s' % run_id
                        ]
                        add_replace_read_groups.extend(option_list)
                        exec_group.add_command(add_replace_read_groups)
