from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class FeatureCounts(AbstractStep):
    '''
    comment here
    '''

    def __init__(self, pipeline):
        super(FeatureCounts, self).__init__(pipeline)

        # in-connections

        # out-connections

        # tools
        self.require_tool('feature_count')

        # required parameters
        self.add_option('a', str, optional=False, description="Give the name \
                        of the annotation file. The program assumes hat the \
                        provided annotation file is in GTF format. Use -F \
                        option to specify other annotation formats.")

        self.add_option('o', str, optional=False, description="Give the name \
                        of the output file. The output file contains the \
                        number of reads assigned to each meta-feature \
                        (or each feature if -f is specified). A meta-feature \
                        is the aggregation of features, grouped by using gene \
                        identifiers. Please refer to the users guide for more \
                        details.")

        # optional parameters
        self.add_option('A', bool, optional=True, default=False,
                        description="Specify the name of a file including \
                        aliases of chromosome names. The file should be a \
                        comma delimited text file that includes two columns. \
                        The first column gives the chromosome names used in \
                        the annotation and the second column gives the \
                        chromosome names used by reads. This file should not \
                        contain header lines. Names included in this file are \
                        case sensitive.")

        self.add_option('F', bool, optional=True, default=False,
                        description="Specify the format of the annotation \
                        file. Acceptable formats include `GTF' and `SAF'. \
                        `GTF' by default. Please refer to the users guide \
                        for SAF annotation format.")

        self.add_option('t', bool, optional=True, default=False,
                        description="Specify the feature type. Only rows \
                        which have the matched feature type in the provided \
                        GTF annotation file will be included for read \
                        counting. `exon' by default.")

        self.add_option('g', bool, optional=True, default=False,
                        description="Specify the attribute type used to group \
                        features (eg. exons) into meta-features (eg. genes), \
                        when GTF annotation is provided. `gene_id' by \
                        default. This attribute type is usually the gene \
                        identifier. This argument is useful for the \
                        meta-feature level summarization.")

        self.add_option('f', bool, optional=True, default=False,
                        description="If specified, read summarization will be \
                        performed at the feature level (eg. exon level). \
                        Otherwise, it is performed at meta-feature level \
                        (eg. gene level).")

        self.add_option('O', bool, optional=True, default=False,
                        description="If specified, reads (or fragments if \
                        -p is specified) will be allowed to be assigned to \
                        more than one matched meta-feature (or feature if \
                        -f is specified).")

        self.add_option('s', int, optional=True, default=0,
                        description="Indicate if strand-specific read \
                        counting should be performed. It has three possible \
                        values:  0 (unstranded), 1 (stranded) and 2 \
                        (reversely stranded). 0 by default.")

        self.add_option('M', bool, optional=True, default=False,
                        description="If specified, multi-mapping \
                        reads/fragments will be counted (ie. a multi-mapping \
                        read will be counted up to N times if it has N \
                        reported mapping locations). The program uses the \
                        `NH' tag to find multi-mapping reads.")

        self.add_option('Q', int, optional=True, default=0,
                        description="The minimum mapping quality score a read \
                        must satisfy in order to be counted. For paired-end \
                        reads, at least one end should satisfy this criteria. \
                        0 by default.")

        self.add_option('T', int, optional=True, default=1,
                        description="Number of the threads. 1 by default.")

        self.add_option('R', bool, optional=True, default=False,
                        description="Output read counting result for each \
                        read/fragment. For each input read file, read \
                        counting results for reads/fragments will be saved \
                        to a tab-delimited file that contains four columns \
                        including read name, status(assigned or the reason \
                        if not assigned), name of target feature/meta-feature \
                        and number of hits if the read/fragment is counted \
                        multiple times. Name of the file is the same as name \
                        of the input read file except a suffix \
                        `.featureCounts' is added.")

        self.add_option('primary', bool, optional=True, default=False,
                        description="If specified, only primary alignments \
                        will be counted. Primary and secondary alignments \
                        are identified using bit 0x100 in the Flag field of \
                        SAM/BAM files. All primary alignments in a dataset \
                        will be counted no matter they are from multi-mapping \
                        reads or not ('-M' is ignored).")

        self.add_option('readExtension5', int, optional=True, default=0,
                        description="Reads are extended upstream by <int> \
                        bases from their 5' end. 0 by default.")

        self.add_option('readExtension3', int, optional=True, default=0,
                        description="Reads are extended upstream by <int> \
                        bases from their 3' end. 0 by default.")

        self.add_option('minReadOverlap', int, optional=True, default=1,
                        description="Specify the minimum number of overlapped \
                        bases required to assign a read to a feature. 1 by \
                        default. Negative values are permitted, indicating a \
                        gap being allowed between a read and a feature.")

        self.add_option('countSplitAlignmentsOnly', bool, optional=True,
                        default=False, description="If specified, only split \
                        alignments (CIGAR strings containing letter `N') will \
                        be counted. All the other alignments will be ignored. \
                        An example of split alignments is the exon-spanning \
                        reads in RNA-seq data.")

        self.add_option('read2pos', int, optional=True, default=0,
                        description="The read is reduced to its 5' most base \
                        or 3' most base. Read summarization is then performed \
                        based on the single base which the read is reduced to.\
                        By default, no read reduction will be performed.")

        self.add_option('ignoreDup', bool, optional=True, default=False,
                        description="If specified, reads that were marked as \
                        duplicates will be ignored. Bit Ox400 in FLAG field \
                        of SAM/BAM file is used for identifying duplicate \
                        reads. In paired end data, the entire read pair will \
                        be ignored if at least one end is found to be a \
                        duplicate read.")

        # Optional paired-end parameters
        self.add_option('p', bool, optional=True, default=False,
                        description="If specified, fragments (or templates) \
                        will be counted instead of reads. This option is \
                        only applicable for paired-end reads. The two reads \
                        from the same fragment must be adjacent to each other \
                        in the provided SAM/BAM file.")

        self.add_option('P', bool, optional=True, default=False,
                        description="If specified, paired-end distance will \
                        be checked when assigning fragments to meta-features \
                        or features. This option is only applicable when -p \
                        is specified. The distance thresholds should be \
                        specified using -d and -D options.")

        self.add_option('d', int, optional=True, default=50,
                        description="Minimum fragment/template length, \
                        50 by default.")

        self.add_option('D', int, optional=True, default=600,
                        description="Maximum fragment/template length, \
                        600 by default.")

        self.add_option('B', bool, optional=True, default=False,
                        description="If specified, only fragments that have \
                        both ends successfully aligned will be considered for \
                        summarization. This option is only applicable for \
                        paired-end reads.")

        self.add_option('C', bool, optional=True, default=False,
                        description="If specified, the chimeric fragments \
                        (those fragments that have their two ends aligned to \
                        different chromosomes) will NOT be included for \
                        summarization. This option is only applicable for \
                        paired-end read data.")

    def runs(self, run_ids_connections_files):
        print('Test')
