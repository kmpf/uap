import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class deepToolsBamCompare(AbstractStep):
    '''
    This tool compares two BAM files based on the number of mapped reads. To
    compare the BAM files, the genome is partitioned into bins of equal size,
    then the number of reads found in each bin is counted per file, and finally
    a summary value is reported. This value can be the ratio of the number of
    reads per bin, the log2 of the ratio, or the difference. This tool can
    normalize the number of reads in each BAM file using the SES method proposed
    by Diaz et al. (2012) “Normalization, bias correction, and peak calling for
    ChIP-seq”. Statistical Applications in Genetics and Molecular Biology, 11(3).
    Normalization based on read counts is also available. The output is either a
    bedgraph or bigWig file containing the bin location and the resulting
    comparison value. By default, if reads are paired, the fragment length
    reported in the BAM file is used. Each mate, however, is treated
    independently to avoid a bias when a mixture of concordant and discordant
    pairs is present. This means that each end will be extended to match the
    fragment length.

    http://deeptools.readthedocs.io/en/latest/content/tools/bamCompare.html

    Usage example::
    
        bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw

    '''

    def __init__(self, pipeline):
        super(deepToolsBamCompare, self).__init__(pipeline)
        
        self.set_cores(10)

        self.add_connection('in/alignments')
        self.add_connection('out/bigwig')

        self.require_tool('bamCompare')

        self.add_option('bed-file', list, optional=True, description='BED file '
                        'that contains all regions that should be considered '
                        'for the coverage analysis. If this option is set '
                        '"multiBamSummary" is executed with "BED-file" '
                        'subcommand, otherwise with "bins" subcommand.')
        self.add_option('samples', dict, optional=False,
                        description='Dictionary with IDs of new runs as keys '
                        'and lists of sample names as values. For each '
                        'sample name a sorted, indexed BAM file is expected '
                        'to be the input from upstream steps.')
        
        # Options for bamCompare
        ## Optional arguments:
        self.add_option('scaleFactorsMethod', str, optional=True,
                        choices=['readCount', 'SES'],
                        description='Method to use to scale the samples. '
                        '(default: readCount)')
        self.add_option('sampleLength', int, optional=True,
                        description="*Only relevant when SES is chosen for the "
                        "scaleFactorsMethod.* To compute the SES, specify the "
                        "length (in bases) of the regions (see "
                        "--numberOfSamples) that will be randomly sampled to "
                        "calculate the scaling factors. If you do not have a "
                        "good sequencing depth for your samples consider "
                        "increasing the sampling regions' size to minimize the "
                        "probability that zero-coverage regions are used. "
                        "(default: 1000)")
        self.add_option('numberOfSamples', int, optional=True,
                        description="*Only relevant when SES is chosen for the "
                        "scaleFactorsMethod.* Number of samplings taken from "
                        "the genome to compute the scaling factors. (default: "
                        "100000.0)")
        self.add_option('scaleFactors', str, optional=True,
                        description='Set this parameter manually to avoid the '
                        'computation of scaleFactors. The format is '
                        'scaleFactor1:scaleFactor2. For example, --scaleFactor '
                        '0.7:1 will cause the first BAM file tobe multiplied '
                        'by 0.7, while not scaling the second BAM file '
                        '(multiplication with 1). (default: None)')
        self.add_option('ratio', str, optional=True,
                        choices=["log2", "ratio", "subtract", "add", "mean",
                                 "reciprocal_ratio", "first,second"],
                        description="The default is to output the log2ratio of "
                        "the two samples. The reciprocal ratio returns the "
                        "negative of the inverse of the ratio if the ratio is "
                        "less than 0. The resulting values are interpreted as "
                        "negative fold changes. *NOTE*: Only with --ratio "
                        "subtract can --normalizeTo1x or --normalizeUsingRPKM "
                        "be used. Instead of performing a computation using "
                        "both files, the scaled signal can alternatively be "
                        "output for the first or second file using the "
                        "'--ratio first' or '--ratio second' (default: log2)")
        self.add_option('pseudocount', float, optional=True,
                        description='small number to avoid x/0. Only useful '
                        'together with --ratio log2 or --ratio ratio . '
                        '(default: 1)')
        self.add_option('binSize', int, optional=True, description='Size of '
                        'the bins, in bases, for the output of the '
                        'bigwig/bedgraph file. (default: 50)')
        # Options for multiBamSummary BED-file subcommand
        ## GTF/BED12 options:
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
        ## Optional arguments:
        self.add_option('blackListFileName', str, optional=True,
                        description="A BED or GTF file containing regions that "
                        "should be excluded from all analyses. Currently this "
                        "works by rejecting genomic chunks that happen to "
                        "overlap an entry. Consequently, for BAM files, if a "
                        "read partially overlaps a blacklisted region or a "
                        "fragment spans over it, then the read/fragment might "
                        "still be considered. Please note that you should "
                        "adjust the effective genome size, if relevant. "
                        "(default: None)")
        self.add_option('region', str, optional=True, description='Region of '
                        'the genome to limit the operation to - this is useful '
                        'when testing parameters to reduce the computing time. '
                        'The format is chr:start:end, for example --region '
                        'chr10 or --region chr10:456700:891000. '
                        '(default: None)')
        ## Output optional options:
        self.add_option('outRawCounts', bool, optional=True,
                        description='Save the counts per region to a '
                        'tab-delimited file. (default: False)')
        ## Read processing options:
        self.add_option('extendReads', int, optional=True,
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
        ## List of options common for bin and BED-file subcommand
        common_options = ['outRawCounts', 'extendReads', 'ignoreDuplicates',
                          'minMappingQuality', 'centerReads', 'samFlagInclude',
                          'samFlagExclude', 'minFragmentLength',
                          'maxFragmentLength', 'region', 'blackListFileName']
        ## List of options for bin subcommand only
        bins_options = ['binSize', 'distanceBetweenBins']
        ## List of options for BED-file subcommand only
        bed_file_options = ['metagene', 'transcriptID',
                            'exonID', 'transcript_id_designator']

        set_options = [option for option in common_options if \
                       self.is_option_set_in_config(option)]

        if self.is_option_set_in_config('bed-file'):
            subcommand = 'BED-file'
            set_options = [option for option in bed_file_options if \
                           self.is_option_set_in_config(option)]
        else:
            subcommand = 'bins'
            set_options = [option for option in bins_options if \
                           self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append( '--%s' % option )
                option_list.append( str(self.get_option(option)) )

        runIds_samples = self.get_option('samples')

        for run_id, samples in runIds_samples.iteritems():
            if not isinstance(run_id, str):
                logger.error("Not a string run ID (%s) for samples (%s)"
                             % (run_id, ", ".join(samples)))
                sys.exit(1)
            if not isinstance(samples, list):
                logger.error("Not a list of samples. Type: %s, Value: %s"
                             % (type(samples), samples))
                sys.exit(1)

            input_paths = list()
            labels = list()
            for sample in samples:
                try:
                    bam_files = run_ids_connections_files[sample]['in/alignments']
                except KeyError:
                    logger.error("No input sample named %s" % sample)
                    sys.exit(1)

                for i in range(len(bam_files)):
                    if not bam_files[i].endswith(".bam"):
                        logger.error("Not a BAM file: %s" % bam_files[i])
                        sys.exit(1)
                    input_paths.append(bam_files[i])
                    if i > 0:
                        labels.append("%s-%s" % (sample, i))
                    else:
                        labels.append(sample)


            with self.declare_run(run_id) as run:
                # Let's compile the command
                with run.new_exec_group() as bam_compare_eg:
                    # 1. multiBamSummary command
                    bam_compare = [
                        self.get_tool('multiBamSummary'),
                        '--bamfile1',
                        '--bamfile2',
                        '--outFileName',
                        run.add_output_file(
                            'bigwig',
                            '%s.bw' % run_id,
                            input_paths
                        )
                    ]
                    bam_compare.extend(input_paths)
                    ## Append output file format
                    bam_compare.append('--outFileFormat')
                    bam_compare.append('bigwig')
                    ## Append list of BED files for BED-file subcommand
                    if subcommand == "BED-file":
                        bam_compare.append('--BED')
                        bam_compare.extend(self.get_option('bed-file'))
                    ## Append list of labels
                    bam_compare.append('--labels')
                    bam_compare.extend(labels)
                    ## Append number of processors
                    bam_compare.extend(['--numberOfProcessors',
                                              str(self.get_cores())])
                    ## Append list of options
                    bam_compare.extend(option_list)

                    bam_compare_eg.add_command(bam_compare)
