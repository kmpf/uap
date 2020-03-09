from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class deepToolsMultiBamSummary(AbstractStep):
    '''
    This step computes the read coverages for genomic regions for every BAM
    input file using multiBamSummary. For downstream analysis such as
    'plotCorrelation' or 'plotPCA' you need to merge the output files.

    The analysis can be performed for the
    entire genome by running the program in 'bins' mode. If you want to count
    the read coverage for specific regions only, use the BED-file mode instead.
    The standard output of multiBamSummary is a compressed numpy array (.npz).
    It can be directly used to calculate and visualize pairwise correlation
    values between the read coverages using the tool 'plotCorrelation'.
    Similarly, plotPCA can be used for principal component analysis of the read
    coverages using the .npz file. Note that using a single bigWig file is only
    recommended if you want to produce a bedGraph file (i.e., with the
    --outRawCounts option; the default output file cannot be used by ANY
    deepTools program if only a single file was supplied!).

    http://deeptools.readthedocs.io/en/latest/content/tools/multiBamSummary.html

    Usage example::

        multiBamSummary [-h] [--version]  ...

    '''

    def __init__(self, pipeline):
        super(deepToolsMultiBamSummary, self).__init__(pipeline)

        self.set_cores(10)

        self.add_connection('in/alignments')
        self.add_connection('out/read-coverage')

        self.require_tool('multiBamSummary')

        self.add_option(
            'bed-file', list, optional=True, description='BED file '
            'that contains all regions that should be considered '
            'for the coverage analysis. If this option is set '
            '"multiBamSummary" is executed with "BED-file" '
            'subcommand, otherwise with "bins" subcommand.')
#        self.add_option('samples', dict, optional=False,
#                        description='Dictionary with IDs of new runs as keys '
#                        'and lists of sample names as values. For each '
#                        'sample name a sorted, indexed BAM file is expected '
#                        'to be the input from upstream steps.')

        # Options for multiBamSummary bins subcommand
        # Optional arguments:
        self.add_option('binSize', int, optional=True, description='Length in '
                        'bases of the window used to sample the genome. '
                        '(default: 10000)')
        self.add_option(
            'distanceBetweenBins',
            int,
            optional=True,
            description='By default, multiBamSummary considers '
            'consecutive bins of the specified --binSize. However, '
            'to reduce the computation time, a larger distance '
            'between bins can by given. Larger distances result in '
            'fewer bins considered. (default: 0)')

        # Options for multiBamSummary BED-file subcommand
        # GTF/BED12 options:
        self.add_option('metagene', bool, optional=True,
                        description="When either a BED12 or GTF file are used "
                        "to provide regions, perform the computation on the "
                        "merged exons, rather than using the genomic interval "
                        "defined by the 5-prime and 3-prime most transcript "
                        "bound (i.e., columns 2 and 3 of a BED file). If a "
                        "BED3 or BED6 file is used as input, then columns 2 "
                        "and 3 are used as an exon. (default: False)")
        self.add_option('transcriptID', str, optional=True,
                        description="When a GTF file is used to provide "
                        "regions, only entries with this value as their "
                        "feature (column 2) will be processed as transcripts. "
                        "(default: transcript)")
        self.add_option('exonID', str, optional=True,
                        description="When a GTF file is used to provide "
                        "regions, only entries with this value as their "
                        "feature (column 2) will be processed as exons. CDS "
                        "would be another common value for this. "
                        "(default: exon)")
        self.add_option('transcript_id_designator', str, optional=True,
                        description="Each region has an ID (e.g., ACTB) "
                        "assigned to it, which for BED files is either column "
                        "4 (if it exists) or the interval bounds. For GTF "
                        "files this is instead stored in the last column as a "
                        "key:value pair (e.g., as 'transcript_id \"ACTB\"', "
                        "for a key of transcript_id and a value of ACTB). In "
                        "some cases it can be convenient to use a different "
                        "identifier. To do so, set this to the desired key. "
                        "(default: transcript_id)")

        # Common options for bins and BED-file subcommand
        # Optional arguments:
        self.add_option(
            'blackListFileName',
            str,
            optional=True,
            description="A BED or GTF file containing regions that "
            "should be excluded from all analyses. Currently this "
            "works by rejecting genomic chunks that happen to "
            "overlap an entry. Consequently, for BAM files, if a "
            "read partially overlaps a blacklisted region or a "
            "fragment spans over it, then the read/fragment might "
            "still be considered. Please note that you should "
            "adjust the effective genome size, if relevant. "
            "(default: None)")
        self.add_option(
            'region', str, optional=True, description='Region of '
            'the genome to limit the operation to - this is useful '
            'when testing parameters to reduce the computing time. '
            'The format is chr:start:end, for example --region '
            'chr10 or --region chr10:456700:891000. '
            '(default: None)')
        # Output optional options:
        self.add_option('outRawCounts', bool, optional=True,
                        description='Save the counts per region to a '
                        'tab-delimited file. (default: False)')
        # Read processing options:
        self.add_option(
            'extendReads',
            bool,
            optional=True,
            description='This parameter allows the extension of '
            'reads to fragment size. If set, each read is '
            'extended, without exception. *NOTE*: This feature is '
            'generally NOT recommended for spliced-read data, such '
            'as RNA-seq, as it would extend reads over skipped '
            'regions. *Single-end*: Requires a user specified '
            'value for the final fragment length. Reads that '
            'already exceed this fragment length will not be '
            'extended. *Paired-end*: Reads with mates are always '
            'extended to match the fragment size defined by the '
            'two read mates. Unmated reads, mate reads that map '
            'too far apart (>4x fragment length) or even map to '
            'different chromosomes are treated like single-end '
            'reads. The input of a fragment length value is '
            'optional. If no value is specified, it is estimated '
            'from the data (mean of the fragment size of all mate '
            'reads). (default: False)')
        self.add_option('ignoreDuplicates', bool, optional=True,
                        description="If set, reads that have the same "
                        "orientation and start position will be considered "
                        "only once. If reads are paired, the mate's position "
                        "also has to coincide to ignore a read. "
                        "(default: False)")
        self.add_option('minMappingQuality', int, optional=True,
                        description="If set, only reads that have a mapping "
                        "quality score of at least this are considered. "
                        "(default: None)")
        self.add_option('centerReads', bool, optional=True,
                        description="By adding this option, reads are centered "
                        "with respect to the fragment length. For paired-end "
                        "data, the read is centered at the fragment length "
                        "defined by the two ends of the fragment. For "
                        "single-end data, the given fragment length is used. "
                        "This option is useful to get a sharper signal around "
                        "enriched regions. (default: False)")
        self.add_option('samFlagInclude', int, optional=True,
                        description="Include reads based on the SAM flag. For "
                        "example, to get only reads that are the first mate, "
                        "use a flag of 64. This is useful to count properly "
                        "paired reads only once, as otherwise the second mate "
                        "will be also considered for the coverage. "
                        "(default: None)")
        self.add_option('samFlagExclude', int, optional=True,
                        description="Exclude reads based on the SAM flag. For "
                        "example, to get only reads that map to the forward "
                        "strand, use --samFlagExclude 16, where 16 is the SAM "
                        "flag for reads that map to the reverse strand. "
                        "(default: None)")
        self.add_option('minFragmentLength', int, optional=True,
                        description="The minimum fragment length needed for "
                        "read/pair inclusion. Note that a value other than 0 "
                        "will exclude all single-end reads. This option is "
                        "primarily useful in ATACseq experiments, for "
                        "filtering mono- or di-nucleosome fragments. "
                        "(default: 0)")
        self.add_option('maxFragmentLength', int, optional=True,
                        description="The maximum fragment length needed for "
                        "read/pair inclusion. A value of 0 disables filtering "
                        "and is needed for including single-end and orphan "
                        "reads. (default: 0)")

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        # List of options common for bin and BED-file subcommand
        common_options = ['outRawCounts', 'extendReads', 'ignoreDuplicates',
                          'minMappingQuality', 'centerReads', 'samFlagInclude',
                          'samFlagExclude', 'minFragmentLength',
                          'maxFragmentLength', 'region']
        # List of options for bin subcommand only
        bins_options = ['binSize', 'distanceBetweenBins']
        # List of options for BED-file subcommand only
        bed_file_options = ['metagene', 'transcriptID', 'exonID',
                            'transcript_id_designator', 'blackListFileName']

        set_options = [option for option in common_options if
                       self.is_option_set_in_config(option)]

        if self.is_option_set_in_config('bed-file'):
            subcommand = 'BED-file'
            set_options.extend([option for option in bed_file_options if
                                self.is_option_set_in_config(option)])
        else:
            subcommand = 'bins'
            set_options.extend([option for option in bins_options if
                                self.is_option_set_in_config(option)])

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():
            # Collect input_paths and labels for multiBamSummary
            input_paths = run_ids_connections_files[run_id]['in/alignments']
            labels = list()
            for f in input_paths:
                if not f.endswith(".bam"):
                    raise StepError(self, "Not a BAM file: %s" % f)
                if len(input_paths) > 1:
                    labels.append("%s-%s" % (run_id, input_paths.index(f)))
                else:
                    labels.append(run_id)

            with self.declare_run(run_id) as run:
                # Let's compile the command
                with run.new_exec_group() as multi_bam_summary_eg:
                    # 1. multiBamSummary command
                    multi_bam_summary = [
                        self.get_tool('multiBamSummary'), subcommand]
                    # Append list of input BAM files
                    multi_bam_summary.append('--bamfiles')
                    multi_bam_summary.extend(input_paths)
                    # Append name of the output file
                    multi_bam_summary.append('--outFileName')
                    multi_bam_summary.append(
                        run.add_output_file(
                            'read-coverage',
                            '%s.npz' % run_id,
                            input_paths
                        )
                    )
                    # Append list of BED files for BED-file subcommand
                    if subcommand == "BED-file":
                        multi_bam_summary.append('--BED')
                        multi_bam_summary.extend(self.get_option('bed-file'))
                    # Append list of labels
                    multi_bam_summary.append('--labels')
                    multi_bam_summary.extend(labels)
                    # Append number of processors
                    multi_bam_summary.extend(['--numberOfProcessors',
                                              str(self.get_cores())])
                    # Append list of options
                    multi_bam_summary.extend(option_list)

                    # Add multiBamSummary to execution group
                    multi_bam_summary_eg.add_command(multi_bam_summary)
