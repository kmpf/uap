from uaperrors import UAPError
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
    by Diaz et al. (2012) "Normalization, bias correction, and peak calling for
    ChIP-seq". Statistical Applications in Genetics and Molecular Biology, 11(3).
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
        self.add_connection('out/ucsc-tracks')

        self.require_tool('bamCompare')

        # Required options:
        self.add_option('samples', list, optional=False,
                        description='List of lists with two elements. '
                        'Each element has to be the name of a run. Each run '
                        'has to provide a SINGLE BAM file. Both BAM files are '
                        'compared using deepTools bamCompare command.')

        # Options for bamCompare
        ## Output options:
        self.add_option('outFileFormat', str, optional=False,
                        choices=["bigwig", "bedgraph"],
                        description='Output file type. Either "bigwig" or '
                        '"bedgraph". (default: "bigwig")')

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
        self.add_option('region', str, optional=True, description='Region of '
                        'the genome to limit the operation to - this is useful '
                        'when testing parameters to reduce the computing time. '
                        'The format is chr:start:end, for example --region '
                        'chr10 or --region chr10:456700:891000. '
                        '(default: None)')
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
        ## Read coverage normalization options:
        self.add_option('normalizeTo1x', int, optional=True,
                        description='Report read coverage normalized to 1x '
                        'sequencing depth (also known as Reads Per Genomic '
                        'Content (RPGC)). Sequencing depth is defined as: '
                        '(total number of mapped reads * fragment length) / '
                        'effective genome size. The scaling factor used is the '
                        'inverse of the sequencing depth computed for the '
                        'sample to match the 1x coverage. To use this option, '
                        'the effective genome size has to be indicated after '
                        'the option. The effective genome size is the portion '
                        'of the genome that is mappable. Large fractions of '
                        'the genome are stretches of NNNN that should be '
                        'discarded. Also, if repetitive regions were not '
                        'included in the mapping of reads, the effective genome '
                        'size needs to be adjusted accordingly. Common values '
                        'are: mm9:2,150,570,000; hg19:2,451,960,000; '
                        'dm3:121,400,000 and ce10:93,260,000. See Table 2 of '
                        'http://www.plosone.org/article/info:doi/10.1371'
                        '/journal.pone.0030377 or http://www.nature.com/nbt'
                        '/journal/v27/n1/fig_tab/nbt.1518_T1.html for several '
                        'effective genome sizes. (default: None)')
        self.add_option('normalizeUsingRPKM', bool, optional=True,
                        description='Use Reads Per Kilobase per Million reads '
                        'to normalize the number of reads per bin. The formula '
                        'is: RPKM (per bin) = number of reads per bin / '
                        '( number of mapped reads (in millions) * bin length '
                        '(kb) ). Each read is considered independently,if you '
                        'want to only count either of the mate pairs in '
                        'paired-end data, use the --samFlag option. (default: '
                        'False)')
        self.add_option('ignoreForNormalization', list, optional=True,
                        description='A list of space-delimited chromosome '
                        'names containing those chromosomes that should be '
                        'excluded for computing the normalization. This is '
                        'useful when considering samples with unequal coverage '
                        'across chromosomes, like male samples. An usage '
                        'examples is --ignoreForNormalization chrX chrM. '
                        '(default: None)')
        self.add_option('skipNonCoveredRegions', bool, optional=True,
                        description='This parameter determines if non-covered '
                        'regions (regions without overlapping reads) in a BAM '
                        'file should be skipped. The default is to treat those '
                        'regions as having a value of zero. The decision to '
                        'skip non-covered regions depends on the interpretation '
                        'of the data. Non-covered regions may represent, for '
                        'example, repetitive regions that should be skipped. '
                        '(default: False)')
        self.add_option('smoothLength', int, optional=True,
                        description='The smooth length defines a window, larger '
                        'than the binSize, to average the number of reads. For '
                        'example, if the --binSize is set to 20 and the '
                        '--smoothLength is set to 60, then, for each bin, the '
                        'average of the bin and its left and right neighbors is '
                        'considered. Any value smaller than --binSize will be '
                        'ignored and no smoothing will be applied. (default: '
                        'None)')
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
        options=['outFileFormat', 'scaleFactorsMethod', 'sampleLength',
                 'numberOfSamples', 'scaleFactors', 'ratio',
                 'pseudocount', 'binSize', 'region', 'blackListFileName',
                 'normalizeTo1x', 'normalizeUsingRPKM',
                 'ignoreForNormalization', 'skipNonCoveredRegions',
                 'smoothLength',
                 'extendReads', 'ignoreDuplicates', 'minMappingQuality',
                 'centerReads', 'samFlagInclude', 'samFlagExclude',
                 'minFragmentLength', 'maxFragmentLength']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]
        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append( '--%s' % option )
                option_list.append( str(self.get_option(option)) )

        # List of sample lists = losl
        losl = self.get_option('samples')

        # Test the user input and connection data for validity
        for samples in losl:
            if len(samples) != 2:
                raise UAPError("Expected exactly two samples. Received %s (%s)"
                             % (len(samples), ", ".join(samples)))
            
            input_paths = list()
            for sample in samples:
                try:
                    files=run_ids_connections_files[sample]['in/alignments']
                except KeyError as e:
                    raise UAPError('No files found for sample %s and connection '
                    '"in/alignments". Please check your configuration.')
                if not len(files) == 1 or not files[0].endswith('.bam'):
                    raise UAPError("Expected exactly one BAM file, got %s"
                                 % ", ".join(files))
                # Add found BAM file to input paths
                input_paths.append(files[0])
            # Assemble new run name from input sample names
            run_id = "%s-%s" % (samples[0], samples[1])
            
            # Start defining the run here:
            with self.declare_run(run_id) as run:
                # Add output file here:
                outfile = str()
                if self.get_option('outFileFormat') == "bigwig":
                    outfile = run.add_output_file(
                        'ucsc-tracks', '%s.bw' % run_id, input_paths)
                elif self.get_option('outFileFormat') == "bedgraph":
                    outfile = run.add_output_file(
                        'ucsc-tracks', '%s.bg' % run_id, input_paths)
                # Let's compile the command
                with run.new_exec_group() as bam_compare_eg:
                    # 1. bamCompare command
                    bam_compare = [
                        self.get_tool('bamCompare'),
                        '--bamfile1', input_paths[0],
                        '--bamfile2', input_paths[1],
                        '--outFileName', outfile
                    ]
                    ## Append number of processors
                    bam_compare.extend(['--numberOfProcessors',
                                              str(self.get_cores())])
                    ## Append list of options
                    bam_compare.extend(option_list)

                    bam_compare_eg.add_command(bam_compare)
