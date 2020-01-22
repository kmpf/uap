import os
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
        self.add_connection(
            'in/alignments',
            #constraints={'min_files_per_run': 1, 'max_files_per_run': 1}
        )
        self.add_connection(
            'in/feature-file',
            #constraints={'total_files': 1}
        )

        # out-connections
        self.add_connection('out/counts')
        self.add_connection('out/summary')
        self.add_connection('out/log_stdout')
        self.add_connection('out/log_stderr')

        # tools
        self.require_tool('feature_counts')

        # required parameters
        self.add_option('a', str, optional=False, description="Give the name \
                        of the annotation file. The program assumes hat the \
                        provided annotation file is in GTF format. Use -F \
                        option to specify other annotation formats.")

        self.add_option('t', str, optional=False, default=None,
                        description="Specify the feature type. Only rows \
                        which have the matched feature type in the provided \
                        GTF annotation file will be included for read \
                        counting. `exon' by default.")

        # optional parameters
        self.add_option('o', str, optional=True, default="counts.txt",
                        description="Give the name \
                        of the output file. The output file contains the \
                        number of reads assigned to each meta-feature \
                        (or each feature if -f is specified). A meta-feature \
                        is the aggregation of features, grouped by using gene \
                        identifiers. Please refer to the users guide for more \
                        details.")

        self.add_option('A', str, optional=True, default=None,
                        description="Specify the name of a file including \
                        aliases of chromosome names. The file should be a \
                        comma delimited text file that includes two columns. \
                        The first column gives the chromosome names used in \
                        the annotation and the second column gives the \
                        chromosome names used by reads. This file should not \
                        contain header lines. Names included in this file are \
                        case sensitive.")

        self.add_option('F', str, optional=True, default=None,
                        description="Specify the format of the annotation \
                        file. Acceptable formats include `GTF' and `SAF'. \
                        `GTF' by default. Please refer to the users guide \
                        for SAF annotation format.")

        self.add_option('g', str, optional=True, default=None,
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

        self.add_option('s', int, optional=False, default=None,
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

        self.add_option('Q', int, optional=True, default=None,
                        description="The minimum mapping quality score a read \
                        must satisfy in order to be counted. For paired-end \
                        reads, at least one end should satisfy this criteria. \
                        0 by default.")

        self.add_option('T', int, optional=True, default=None,
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

        self.add_option('readExtension5', int, optional=True, default=None,
                        description="Reads are extended upstream by <int> \
                        bases from their 5' end. 0 by default.")

        self.add_option('readExtension3', int, optional=True, default=None,
                        description="Reads are extended upstream by <int> \
                        bases from their 3' end. 0 by default.")

        self.add_option('minReadOverlap', int, optional=True, default=None,
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

        self.add_option('read2pos', int, optional=True, default=None,
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

        self.add_option('d', int, optional=True, default=None,
                        description="Minimum fragment/template length, \
                        50 by default.")

        self.add_option('D', int, optional=True, default=None,
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

        self.set_cores(4)
        self.add_option('cores', int, default=3)

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                alignments = run_ids_connections_files[run_id]['in/alignments']
                input_paths = alignments
                feature_path = str

                # try to get feature file from in_connection or option
                try:
                    feature_path = run_ids_connections_files[run_id]['in/feature-file'][0]
                except KeyError:
                    if self.is_option_set_in_config('a'):
                        feature_path = self.get_option('a')
                    else:
                        logger.error(
                            "No feature file could be found for '%s'" % run_id)
                        exit(1)
                    if not os.path.isfile(feature_path):
                        logger.error("Feature file '%s' is not a file."
                                     % feature_path)
                        exit(1)

                fc_exec_group = run.new_exec_group()
                fc = [self.get_tool('feature_counts')]

                if self.is_option_set_in_config('A'):
                    fc.extend(['-A', str(self.get_option('A'))])

                if self.is_option_set_in_config('F'):
                    fc.extend(['-F', str(self.get_option('F'))])

                if self.is_option_set_in_config('t'):
                    fc.extend(['-t', str(self.get_option('t'))])

                if self.is_option_set_in_config('g'):
                    fc.extend(['-g', str(self.get_option('g'))])

                if self.get_option('f') is True:
                    fc.extend(['-f'])

                if self.get_option('O') is True:
                    fc.extend(['-O'])

                if self.is_option_set_in_config('s'):
                    fc.extend(['-s', str(self.get_option('s'))])

                if self.get_option('M') is True:
                    fc.extend(['-M'])

                if self.is_option_set_in_config('Q'):
                    fc.extend(['-Q', str(self.get_option('Q'))])

                if self.is_option_set_in_config('T'):
                    fc.extend(['-T', str(self.get_option('T'))])

                if self.get_option('R') is True:
                    fc.extend(['-R'])

                if self.get_option('primary') is True:
                    fc.extend(['--primary'])

                if self.is_option_set_in_config('readExtension5'):
                    fc.extend(['--readExtension5',
                               str(self.get_option('readExtension5'))])

                if self.is_option_set_in_config('readExtension3'):
                    fc.extend(['--readExtension3',
                               str(self.get_option('readExtension3'))])

                if self.is_option_set_in_config('minReadOverlap'):
                    fc.extend(['--minReadOverlap',
                               str(self.get_option('minReadOverlap'))])

                if self.get_option('countSplitAlignmentsOnly') is True:
                    fc.extend(['--countSplitAlignmentsOnly'])

                if self.is_option_set_in_config('minReadOverlap'):
                    fc.extend(['--minReadOverlap',
                               str(self.get_option('minReadOverlap'))])

                if self.get_option('ignoreDup') is True:
                    fc.extend(['--ignoreDup'])

                if self.get_option('p') is True:
                    fc.extend(['-p'])

                if self.get_option('P') is True:
                    fc.extend(['-P'])

                if self.is_option_set_in_config('d'):
                    fc.extend(['-d', str(self.get_option('d'))])

                if self.is_option_set_in_config('D'):
                    fc.extend(['-D', str(self.get_option('D'))])

                if self.get_option('B') is True:
                    fc.extend(['-B'])

                if self.get_option('C') is True:
                    fc.extend(['-C'])

                fc.extend(['-a', os.path.abspath(feature_path)])

                basename = run_id + '.' +  self.get_option('o')
                fc.extend(['-o', basename])

                fc.extend(input_paths)

                stderr_file = "%s-featureCounts_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, input_paths)
                stdout_file = "%s-featureCounts_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, input_paths)

                run.add_output_file('counts',
                                    run_id + '.' +  self.get_option('o'),
                                    input_paths)
                run.add_output_file('summary',
                                    run_id + '.' + self.get_option('o') + '.summary',
                                    input_paths)
                fc_exec_group = run.new_exec_group()
                fc_exec_group.add_command(fc, stdout_path=log_stdout,
                                          stderr_path=log_stderr)
