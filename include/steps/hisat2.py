from uaperrors import UAPError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Hisat2(AbstractStep):
    '''
    HISAT2 is a fast and sensitive alignment program for mapping
    next-generation sequencing reads (both DNA and RNA) to a population of
    human genomes (as well as to a single reference genome).

    https://ccb.jhu.edu/software/hisat2/index.shtml
    must be version 2.1 or higher
    metrics and summary file are automatically produced
    '''

    def __init__(self, pipeline):
        super(Hisat2, self).__init__(pipeline)
        self.set_cores(12)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read', optional=True)
        self.add_connection('out/alignments')
        self.add_connection('out/log_stderr')
        self.add_connection('out/metrics')
        self.add_connection('out/summary')
        self.add_connection('out/unaligned', optional=True, format='fastq.gz',
                description='Unpaired reads that didn\'t align.')
        self.add_connection('out/aligned', optional=True, format='fastq.gz',
                description='Unpaired reads that aligned.')

        self.require_tool('pigz')
        self.require_tool('hisat2')

        self.add_option('index', str, optional=False,
                        description="Path to hisat2 index (not containing \
                        file suffixes).")

        self.add_option('cores', int, default=12)

        # Input:
        self.add_option('q', bool, default=False, optional=True,
                        description="query input files are FASTQ .fq/.fastq \
                        (default)")

        self.add_option('qseq', bool, default=False, optional=True,
                        description="query input files are in Illumina's \
                        qseq format")

        self.add_option('f', bool, default=False, optional=True,
                        description="query input files are (multi-)FASTA \
                        .fa/.mfa")

        self.add_option('r', bool, default=False, optional=True,
                        description="query input files are \
                        raw one-sequence-per-line")

        self.add_option('c', bool, default=False, optional=True,
                        description="<m1>, <m2>, <r> are sequences \
                        themselves, not files")

        self.add_option('skip', int, optional=True,
                        description="skip the first <int> reads/pairs \
                        in the input (none)")
        self.add_option('upto', int, optional=True,
                        description="stop after first <int> reads/pairs \
                        (no limit)")
        self.add_option('trim5', int, optional=True,
                        description="trim <int> bases from 5'/left end \
                        of reads (0)")
        self.add_option('trim3', int, optional=True,
                        description="trim <int> bases from 3'/right end \
                        of reads (0)")
        self.add_option('phred33', bool, optional=True, default=False,
                        description="qualities are Phred+33 (default)")
        self.add_option('phred64', bool, optional=True, default=False,
                        description="qualities are Phred+64")
        self.add_option('int-quals', bool, optional=True, default=False,
                        description="qualities encoded as space-delimited \
                        integers")

        # Presets:?

        # Alignment:

        # N, L, i?

        self.add_option('ignore-quals', bool, default=None, optional=True,
                        description="treat all quality values as 30 on Phred \
                        scale (off)")

        self.add_option('n-ceil', str, default=None, optional=True,
                        description="func for max # non-A/C/G/Ts permitted in \
                        aln (L,0,0.15)")

        # dpad, gbar?

        self.add_option('nofw', bool, default=None, optional=True,
                        description="do not align forward (original) version \
                        of read (off)")

        self.add_option('norc', bool, default=None, optional=True,
                        description="do not align reverse-complement version \
                        of read (off)")

        # Spliced Alignment:
        self.add_option('pen-cansplice', str, default=None, optional=True,
                        description="penalty for a canonical splice site (0)")

        self.add_option('pen-noncansplice', str, default=None, optional=True,
                        description="penalty for a non-canonical splice site \
                        (12)")

        self.add_option('pen-canintronlen', str, default=None, optional=True,
                        description="penalty for long introns (G,-8,1) with \
                        canonical splice sites")

        self.add_option('pen-noncanintronlen', str, default=None,
                        optional=True, description="penalty for long introns \
                        (G,-8,1) with noncanonical splice sites")

        self.add_option('min-intronlen', str, default=None, optional=True,
                        description="minimum intron length (20)")

        self.add_option('max-intronlen', str, default=None, optional=True,
                        description="maximum intron length (500000)")

        self.add_option('known-splicesite-infile', str, default=None,
                        optional=True, description="provide a list of known \
                        splice sites")

        self.add_option('novel-splicesite-outfile', str, default=None,
                        optional=True, description="report a list of splice \
                        sites")

        self.add_option('novel-splicesite-infile', str, default=None,
                        optional=True, description="provide a list of novel \
                        splice sites")

        self.add_option('no-temp-splicesite', bool, default=False,
                        optional=True, description="disable the use of splice \
                        sites found")

        self.add_option('no-spliced-alignment', bool, default=False,
                        optional=True, description="disable spliced alignment")

        # note to self just allow R F U and if paired in extend accordingly
        # otherwise mixed single paired will not work
        # Truseq is RF (R)
        self.add_option('rna-strandness', str, choices=["R", "F", "U"],
                        default=None, optional=False,
                        description="Specify strand-specific information \
                        (unstranded); paired and are extended F->FR, R->RF")

        self.add_option('tmo', bool, default=False, optional=True,
                        description="Reports only those alignments within \
                        known transcriptome")

        self.add_option('dta', bool, default=False, optional=True,
                        description="Reports alignments tailored for \
                        transcript assemblers")

        # dta-cufflinks?

        # Scoring:

        self.add_option('ma', str, default=None, optional=True,
                        description="match bonus (0 for --end-to-end, 2 for \
                        --local)")

        self.add_option('mp', str, default=None, optional=True,
                        description="max and min penalties for mismatch; \
                        lower qual = lower penalty <2,6>")

        self.add_option('sp', str, default=None, optional=True,
                        description="max and min penalties for soft-clipping; \
                        lower qual = lower penalty <1,2>")

        self.add_option('no-softclip', bool, optional=True,
                        description='no soft-clipping')

        self.add_option('np', str, default=None, optional=True,
                        description="penalty for non-A/C/G/Ts in read/ref (1)")

        self.add_option('rdg', str, default=None, optional=True,
                        description="read gap open, extend penalties (5,3)")

        self.add_option('rfg', str, default=None, optional=True,
                        description="reference gap open, extend penalties \
                        (5,3)")

        self.add_option('score-min', str, default=None, optional=True,
                        description="min acceptable alignment score w/r/t \
                        read length (G,20,8 for local, L,-0.6,-0.6 for \
                        end-to-end)")

        # Reporting:

        self.add_option('k', int, default=None, optional=True,
                        description="report up to <int> alns per read; \
                        MAPQ not meaningful")

        # a?

        # Effort:

        # D, R?

        ###############
        # Paired-end: #
        ###############

        # I, X?

        self.add_option('minins', int, optional=True,
                        description='minimum fragment length (0), only valid with '
                        '--no-spliced-alignment')

        self.add_option('maxins', int, optional=True,
                        description='maximum fragment length (500), only valid with '
                        '--no-spliced-alignment')

        self.add_option('library_type', str, optional=False,
                        choices=['fr', 'rf', 'ff'], default='fr',
                        description='-1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)')

        self.add_option('no-mixed', bool, default=None, optional=True,
                        description="suppress unpaired alignments for paired \
                        reads")

        self.add_option('no-discordant', bool, default=None, optional=True,
                        description="suppress discordant alignments for \
                        paired reads")

        # no-dovetail, no-contain, no-overlap?

        ################
        # Performance: #
        ################

        self.add_option('offrate', int, optional=True,
                        description='override offrate of index; must be >= index\'s offrate')
        self.add_option('reorder', bool, optional=True,
                        description='force SAM output order to match order of input reads')
        self.add_option('mm', bool, optional=True,
                        description='use memory-mapped I/O for index; many \'hisat2\'s can share')

        ###########
        # Output: #
        ###########

        # t, un, al, un-conc, al-conc?

        self.add_option('un-gz', bool, default=False, optional=True,
                        description='write unpaired reads that didn\'t align \
                        to gzip compress output connection "out/unaligned"')

        self.add_option('al-gz', bool, default=False, optional=True,
                        description='write unpaired reads that aligned \
                        to gzip compress output connection "out/aligned"')

        self.add_option('quiet', bool, default=False, optional=True,
                        description="print nothing to stderr except serious \
                        errors")

        # met-file, met-stderr, met?

        self.add_option('new-summary', bool, default=False, optional=True,
                        description="print alignment summary in a new style, \
                        which is more machine-friendly")

        self.add_option('no-head', bool, default=False, optional=True,
                        description="supppress header lines, i.e. lines \
                        starting with @")

        self.add_option('no-sq', bool, default=False, optional=True,
                        description="supppress @SQ header lines")

        self.add_option('rg-id', str, default=None, optional=True,
                        description="puts sample name in rg")

        self.add_option('rg', str, default=None, optional=True,
                        description="add <text> ('lab:value') to @RG line of \
                        SAM header. (Note: @RG line only printed when --rg-id \
                        is set.)")

        self.add_option('omit-sec-seq', bool, default=False, optional=True,
                        description="put '*' in SEQ and QUAL fields for \
                        secondary alignments")

        # notice: params add-chrname and remove-chrname
        # available from version 2.0.4
        self.add_option('add-chrname', bool, default=False, optional=True,
                        description="Add 'chr' to reference names in \
                        alignment (e.g., 18 to chr18)")

        self.add_option('remove-chrname', bool, default=False, optional=True,
                        description="Remove 'chr' from reference names in \
                        alignment (e.g., chr18 to 18)")

        ################
        # Performance: #
        ################

        ##########
        # Other: #
        ##########

        self.add_option('qc-filter', bool, default=False, optional=True,
                        description="filter out reads that are bad according \
                        to QSEQ filter")

        self.add_option('seed', int, optional=True,
                        description='seed for random number generator (0)')

        self.add_option('non-deterministic', bool, default=False, optional=True,
                        description="seed rand. gen. arbitrarily instead of \
                        using read attributes")

    def runs(self, cc):
        flags = ['q', 'qseq', 'skip', 'f', 'c', 'ignore-quals', 'nofw', 'dta',
                 'norc', 'no-mixed',  'no-discordant', 'quiet', 'qc-filter',
                 'non-deterministic', 'no-temp-splicesite', 'no-softclip',
                 'no-spliced-alignment', 'tmo', 'no-head', 'no-sq',
                 'omit-sec-seq', 'remove-chrname', 'add-chrname', 'new-summary']

        strflags = ['n-ceil', 'ma', 'mp', 'sp', 'np', 'rdg', 'score-min', 'k',
                    'skip', 'rfg', 'rg', 'pen-cansplice', 'pen-noncansplice',
                    'pen-canintronlen', 'pen-noncanintronlen', 'min-intronlen',
                    'max-intronlen', 'known-splicesite-infile',
                    'minins', 'maxins', 'seed', 'trim5', 'trim3',
                    'novel-splicesite-outfile', 'novel-splicesite-infile']

        self.set_cores(self.get_option('cores'))

        # Check if option values are valid
        if not os.path.exists(self.get_option('index') + '.1.ht2'):
            raise UAPError("Could not find index file: %s.*" %
                         self.get_option('index'))

        paired_end = cc.connection_exists('in/second_read')

        if not cc.all_runs_have_connection('in/first_read'):
            read_name = '' if paired_end else ' first'
            run_ids = list(cc.get_runs_without_any('in/first_read'))
            if len(run_ids)>5:
                run_ids = run_ids[0:5] + ['...']
            raise UAPError('[Hisat2] No%s read passed by runs '
                           '%s.' % (read_name, list(run_ids)))

        if paired_end and not cc.all_runs_have_connection('in/second_read'):
            run_ids = list(cc.get_runs_without_any('in/second_read'))
            if len(run_ids)>5:
                run_ids = run_ids[0:5] + ['...']
            raise UAPError('[Hisat2] No second read passed by runs '
                           '%s.' % run_ids)

        for run_id in cc.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = cc[run_id]['in/first_read'][0]
                input_paths = [fr_input]
                is_paired_end = False
                if paired_end:
                    sr_input = cc[run_id]['in/second_read'][0]
                    input_paths.append(sr_input)

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as hisat2_pipe:
                        # Assemble hisat2 command
                        hisat2 = [self.get_tool('hisat2')]

                        for flag in flags:
                            if self.is_option_set_in_config(flag):
                                if self.get_option(flag) is True:
                                    if flag in ['q', 'f', 'r', 'c']:
                                        hisat2.extend(['-' + flag])
                                    else:
                                        hisat2.extend(['--' + flag])

                        lt = self.get_option('library_type')
                        hisat2.append('--%s' % lt)

                        for flag in strflags:
                            if self.is_option_set_in_config(flag):
                                hisat2.extend(['--' + flag,
                                              str(self.get_option(flag))])

                        # Leave 2 cores available for pigz compressing the output.
                        hisat2.extend(['-p', str(self.get_option('cores') - 2),
                                       '-x', os.path.abspath(self.get_option('index'))])

                        if paired_end:
                            if self.get_option('rna-strandness') == 'F':
                                hisat2.extend(['--rna-strandness', 'FR'])
                            elif self.get_option('rna-strandness') == 'R':
                                hisat2.extend(['--rna-strandness', 'RF'])

                            hisat2.extend([
                                '-1', fr_input,
                                '-2', sr_input])
                        else:
                            hisat2.extend(['-U', fr_input])
                            if not self.get_option('rna-strandness') == 'U':
                                hisat2.extend([
                                    '--rna-strandness',
                                    self.get_option('rna-strandness')])

                        log_stderr = run.add_output_file(
                            'log_stderr',
                            '%s-hisat2-log_stderr.txt' % run_id,
                            input_paths)

                        summary = run.add_output_file(
                            'summary',
                            '%s-hisat2-summary.txt' % run_id,
                            input_paths)
                        hisat2.extend(['--summary-file', summary])

                        metrics = run.add_output_file(
                            'metrics',
                            '%s-hisat2-metrics.txt' % run_id,
                            input_paths)
                        hisat2.extend(['--met-file', metrics])

                        if self.get_option('un-gz') is True:
                            unaligned = run.add_output_file(
                                'unaligned',
                                '%s-hisat2-unaligned.fastq.gz' % run_id,
                                input_paths)
                            hisat2.extend(['--un-gz', unaligned])

                        if self.get_option('al-gz') is True:
                            unaligned = run.add_output_file(
                                'aligned',
                                '%s-hisat2-aligned.fastq.gz' % run_id,
                                input_paths)
                            hisat2.extend(['--al-gz', aligned])

                        hisat2_pipe.add_command(hisat2, stderr_path=log_stderr)
                        res = run.add_output_file(
                            'alignments',
                            '%s-hisat2-results.sam.gz' % run_id,
                            input_paths)

                        # Compress hisat2 output
                        pigz = [self.get_tool('pigz'),
                                '--stdout']
                        hisat2_pipe.add_command(pigz, stdout_path=res)
