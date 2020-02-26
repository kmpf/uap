from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class deepToolsPlotFingerprint(AbstractStep):
    '''
    This tool samples indexed BAM files and plots a profile of cumulative
    read coverages for each. All reads overlapping a window (bin) of the
    specified length are counted; these counts are sorted and the cumulative
    sum is finally plotted.

    Usage example::

       plotFingerprint -b treatment.bam control.bam -plot fingerprint.png

    '''
    def __init__(self, pipeline):
        super(deepToolsPlotFingerprint, self).__init__(pipeline)

        self.set_cores(10)

        self.add_connection('in/alignments')
        self.add_connection('out/plots')
        self.add_connection('out/counts')

        self.require_tool('plotFingerprint')

        # Required options:
        self.add_option('samples', list, optional=False,
                        description='List of lists with run names. '
                        'Each element has to be the name of a run. Each run '
                        'has to provide a SINGLE BAM file. All BAM files are '
                        'plotted and counted using deepTools plotFingerprint '
                        'command.')

        # Output options:
        self.add_option('plotFileFormat', str, optional=False,
                        choices=["png", "eps", "pdf", "svg"],
                        description="File ending of the output figure. It will "
                        "be used to determine the image format.")

        # Read processing options:
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

        # Optional arguments:
        self.add_option('binSize', int, optional=True,
                        description="Window size in base pairs to sample the "
                        "genome. (default: 500)")
        self.add_option("numberOfSamples", int,  optional=True,
                        description="Number of bins that sampled from the "
                        "genome, for which the overlapping number of reads is "
                        "computed. (default: 500000.0)")
        #--plotTitle PLOTTITLE, -T PLOTTITLE
        #                Title of the plot, to be printed on top of the generated
        #                image. Leave blank for no title. (default: )
        self.add_option("skipZeros", bool, optional=True,
                        description="If set, then regions with zero overlapping "
                        "reads for *all* given BAM files are ignored. This will "
                        "result in a reduced number of read counts than that "
                        "specified in --numberOfSamples "
                        "(default: False)")
        self.add_option("outQualityMetrics", str, optional=True,
                        description="Quality metrics can optionally be output "
                        "to this file. The file will have one row per input BAM "
                        "file and columns containing a number of metrics. "
                        "Please see the online documentation for a longer "
                        "explanation: http://deeptools.readthedocs.io/en/latest"
                        "/content/feature/plotFingerprint_QC_metrics.html. "
                        "(default: None)")
        self.add_option("JSDsample", str, optional=True,
                        description='Reference sample against which to compute '
                        'the Jensen-Shannon distance and the CHANCE statistics. '
                        'If this is not specified, then these will not be '
                        'calculated. If --outQualityMetrics is not specified '
                        'then this will be ignored. The Jensen-Shannon '
                        'implementation is based on code from Sitanshu Gakkhar '
                        'at BCGSC. The CHANCE implementation is based on code '
                        'from Matthias Haimel. (default: None)')
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

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        ## List of options common for bin and BED-file subcommand
        options=['extendReads', 'ignoreDuplicates', 'minMappingQuality',
                 'centerReads', 'samFlagInclude', 'samFlagExclude',
                 'minFragmentLength', 'maxFragmentLength', 'binSize',
                 'numberOfSamples', 'skipZeros', 'outQualityMetrics',
                 'JSDsample', 'region', 'blackListFileName']

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

        losl = self.get_option('samples')

        for samples in losl:
            run_id = str()
            input_paths = list()
            labels = list()
            for sample in samples:
                bamfiles = run_ids_connections_files[sample]['in/alignments']
                if not len(bamfiles) == 1:
                    raise StepError(self, "Expected a single file for input run %s "
                                 "received %s" % (sample, bamfiles))
                if not bamfiles[0].endswith(".bam"):
                    raise StepError(self, "Not a BAM file: %s" % bamfiles[0])

                input_paths.append(bamfiles[0])
                labels.append(sample)

            run_id = "-".join(labels)
            with self.declare_run(run_id) as run:
                # Let's compile the command
                with run.new_exec_group() as plot_fingerprint_eg:
                    plotfile = run.add_output_file(
                        'plots',
                        '%s.%s' % (run_id,
                                   self.get_option('plotFileFormat')),
                        input_paths)
                    countsfile = run.add_output_file(
                        'counts',
                        '%s.counts' % run_id,
                        input_paths)

                    plot_fingerprint = [
                        self.get_tool("plotFingerprint"),
                        '--plotFile', plotfile,
                        '--outRawCounts', countsfile,
                        '--numberOfProcessors', str(self.get_cores()),
                        '--bamfiles']
                    plot_fingerprint.extend(input_paths)
                    plot_fingerprint.append('--labels')
                    plot_fingerprint.extend(labels)

                    plot_fingerprint.extend(option_list)

                    plot_fingerprint_eg.add_command(plot_fingerprint)
