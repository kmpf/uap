from uaperrors import UAPError
import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PicardCollectRnaSeqMetrics(AbstractStep):
    '''
    Produces RNA alignment metrics for a SAM or BAM file.

    This tool takes a SAM/BAM file containing the aligned reads from an RNAseq experiment and produces metrics describing the distribution of the bases within the transcripts. It calculates the total numbers and the fractions of nucleotides within specific genomic regions including untranslated regions (UTRs), introns, intergenic sequences (between discrete genes), and peptide-coding sequences (exons). This tool also determines the numbers of bases that pass quality filters that are specific to Illumina data (PF_BASES). For more information please see the corresponding GATK Dictionary entry.

    Other metrics include the median coverage (depth), the ratios of 5 prime /3 prime-biases, and the numbers of reads with the correct/incorrect strand designation. The 5 prime /3 prime-bias results from errors introduced by reverse transcriptase enzymes during library construction, ultimately leading to the over-representation of either the 5 prime or 3 prime ends of transcripts. Please see the CollectRnaSeqMetrics definitions for details on how these biases are calculated.

    The sequence input must be a valid SAM/BAM file containing RNAseq data aligned by an RNAseq-aware genome aligner such a STAR or TopHat. The tool also requires a REF_FLAT file, a tab-delimited file containing information about the location of RNA transcripts, exon start and stop sites, etc. For more information on the REF_FLAT format, see the following description. Build-specific REF_FLAT files can be obtained here.
    Usage example:

    java -jar picard.jar CollectRnaSeqMetrics \
    	   I=input.bam \
           O=output.RNA_Metrics \
	   REF_FLAT=ref_flat.txt \
	   STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
	   RIBOSOMAL_INTERVALS=ribosomal.interval_list

       https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics

    '''

    def __init__(self, pipeline):
        super(PicardCollectRnaSeqMetrics, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/alignments')
	self.add_connection('in/refFlat')
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

        # [Picard CollectRnaSeqMetrics Options:]
        self.add_option('REF_FLAT', str, optional = False)
        self.add_option('RIBOSOMAL_INTERVALS', str, optional = True)
        self.add_option('STRAND_SPECIFICITY', str, choices=['NONE', 'FIRST_READ_TRANSCRIPTION_STRAND', 'SECOND_READ_TRANSCRIPTION_STRAND'],  optional = False)
        self.add_option('MINIMUM_LENGTH', int, optional = True)
        self.add_option('IGNORE_SEQUENCE', str, optional = True) # not implemented 

        self.add_option('ASSUME_SORTED', bool, optional = True)
        self.add_option('RRNA_FRAGMENT_PERCENTAGE', int, optional = True)
        self.add_option('METRIC_ACCUMULATION_LEVEL', str, choices=['ALL_READS', 'SAMPLE', 'LIBRARY', 'READ_GROUP'], optional = True)
        self.add_option('STOP_AFTER', str, optional = True)


    def runs(self, run_ids_connections_files):

        options = [
            # Standard Picard Options:
            'TMP_DIR', 'VERBOSITY', 'QUIET', 'VALIDATION_STRINGENCY',
            'COMPRESSION_LEVEL', 'MAX_RECORDS_IN_RAM', 'CREATE_INDEX',
            'CREATE_MD5_FILE', 'REFERENCE_SEQUENCE', 'GA4GH_CLIENT_SECRETS',
            # Picard CollectRnaSeqMetrics Options:
            'REF_FLAT', 'RIBOSOMAL_INTERVALS', 'STRAND_SPECIFICITY', 'MINIMUM_LENGTH',
            'IGNORE_SEQUENCE', 'ASSUME_SORTED', 'RRNA_FRAGMENT_PERCENTAGE',
            'METRIC_ACCUMULATION_LEVEL', 'STOP_AFTER']
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
                        metrics = run.add_output_file(
                            "metrics", '%s-picard-rna-seq-metrics.txt' % run_id,
                            input_paths)

                        collectM = [
                            self.get_tool('picard-tools'), 'CollectRnaSeqMetrics',
                            'INPUT=%s' % input_paths[0],
                            'OUTPUT=%s' % metrics,
                        ]
                        collectM.extend(option_list)
                        exec_group.add_command(collectM)
