import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class deepToolsBamPEFragmentSize(AbstractStep):
    '''
    bamPEFragmentSize tool calculates the fragment sizes for read pairs given a
    BAM file from paired-end sequencing.Several regions are sampled depending on
    the size of the genome and number of processors to estimate thesummary
    statistics on the fragment lengths. Properly paired reads are preferred for
    computation, i.e., it will only use discordant pairs if no concordant
    alignments overlap with a given region. The default setting simply prints
    the summary statistics to the screen.

    http://deeptools.readthedocs.io/en/latest/content/tools/bamPEFragmentSize.html

    Usage example::
    
        bamPEFragmentSize [-h] [--bamfiles bam files [bam files ...]]
                          [--histogram FILE] [--numberOfProcessors INT]
                          [--samplesLabel SAMPLESLABEL [SAMPLESLABEL ...]]
                          [--plotTitle PLOTTITLE]
                          [--maxFragmentLength MAXFRAGMENTLENGTH] [--logScale]
                          [--binSize INT] [--distanceBetweenBins INT]
                          [--blackListFileName BED file] [--verbose] [--version]
 
    '''

    def __init__(self, pipeline):
        super(deepToolsBamPEFragmentSize, self).__init__(pipeline)
        
        self.set_cores(10)

        self.add_connection('in/alignments')
        self.add_connection('out/fragment_size_stats')
        self.add_connection('out/fragment_size_plots')

        self.require_tool('bamPEFragmentSize')

        self.add_option('samples', dict, optional=True,
                        description='Dictionary with IDs of new runs as keys '
                        'and lists of sample names as values. For each sample '
                        'name a BAM file is expected to be the input from '
                        'upstream steps. If not provided this step '
                        'calculates summary statistics for each input file.')
        self.add_option('histogram', bool, optional=True,
                        description='If set saves a .png file with a histogram '
                        'of fragment length distribution for each run.')
        self.add_option('maxFragmentLength', int, optional=True,
                        description='The maximum fragment length in the '
                        'histogram. A value of 0 (the default) indicates to '
                        'use twice the mean fragment length')
        self.add_option('logScale', bool, optional=True,
                        description='Plot on the log scale')
        self.add_option('binSize', int, optional=True,
                        description='Length in bases of the window used to '
                        'sample the genome. (default 1000)')
        self.add_option('distanceBetweenBins', int, optional=True,
                        description='To reduce the computation time, not every '
                        'possible genomic bin is sampled. This option allows '
                        'you to set the distance between bins actually sampled '
                        'from. Larger numbers are sufficient for high coverage '
                        'samples, while smaller values are useful for lower '
                        'coverage samples. Note that if you specify a value '
                        'that results in too few (<1000) reads sampled, the '
                        'value will be decreased. (default 1000000)')
        self.add_option('blackListFileName', str, optional=True,
                        description="A BED file containing regions that "
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
        options = ['samplesLabel', 'maxFragmentLength',
                   'binSize', 'distanceBetweenBins',
                   'blackListFileName']
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

        def declare_bamPEFragmentSize(run_id, input_paths, labels):
            with self.declare_run(run_id) as run:
                # Let's compile the command
                with run.new_exec_group() as bamPEFragmentSize_eg:
                    # 1. bamPEFragmentSize command
                    bamPEFragmentSize = [
                        self.get_tool('bamPEFragmentSize'),
                        '--numberOfProcessors', self.get_cores(),
                        '--bamfiles']
                    bamPEFragmentSize.extend(input_paths)

                    # Set options for plot creation
                    if self.is_option_set_in_config('histogram'):
                        bamPEFragmentSize.append('--histogram')
                        bamPEFragmentSize.append(
                            run.add_output_file(
                                'fragment_size_plots',
                                '%s.png' % run_id,
                                input_paths))
                        bamPEFragmentSize.append('--plotTitle')
                        bamPEFragmentSize.append(run_id)
                        bamPEFragmentSize.append('--samplesLabel')
                        bamPEFragmentSize.extend(labels)
                        if self.is_option_set_in_config('logScale'):
                            bamPEFragmentSize.append('--logScale')
                        
                    ## Append list of options
                    bamPEFragmentSize.extend(option_list)

                    bamPEFragmentSize_eg.add_command(
                        bamPEFragmentSize,
                        stdout_path = run.add_output_file(
                            'fragment_size_stats',
                            '%s-PEFragmentSize.stats' % run_id,
                            input_paths))


                
        run_id = str()
        input_paths = list()
        labels = list()
        if self.is_option_set_in_config('samples'):
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
                # Start declaring the command
                declare_bamPEFragmentSize(run_id, input_paths, labels)

        else:
            for run_id in run_ids_connections_files.keys():
                try:
                    input_paths = run_ids_connections_files[run_id]['in/alignments']
                except KeyError:
                    logger.error("No input sample named %s" % sample)
                    sys.exit(1)
                for f in input_paths:
                    label = run_id
                    if len(input_paths) > 1:
                        label = "%s-%s" % (run_id, input_paths.index(f))
                    declare_bamPEFragmentSize(run_id, f, [label])
