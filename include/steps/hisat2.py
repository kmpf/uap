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
    '''

    def __init__(self, pipeline):
        super(Hisat2, self).__init__(pipeline)
        self.set_cores(12)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/log_stderr')
        self.add_connection('out/metrics')
        self.add_connection('out/unaligned')

        self.require_tool('pigz')
        self.require_tool('hisat2')

        self.add_option('index', str, optional=False,
                        description="Path to hisat2 index (not containing \
                        file suffixes).")

        self.add_option('cores', int, default=12)

        # TODO: remove-chrname, add-chrname not in help list

        # Input:
        self.add_option('q', bool, default=None, optional=True,
                        description="query input files are FASTQ .fq/.fastq \
                        (default)")

        self.add_option('qseq', bool, default=None, optional=True,
                        description="query input files are in Illumina's \
                        qseq format")

        self.add_option('f', bool, default=None, optional=True,
                        description="query input files are (multi-)FASTA \
                        .fa/.mfa")

        # r?

        self.add_option('c', bool, default=None, optional=True,
                        description="<m1>, <m2>, <r> are sequences \
                        themselves, not files")

        # s, u, 5, 3, phred33, phred64, int-quals

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

        self.add_option('no-temp-splicesite', bool, default=None,
                        optional=True, description="disable the use of splice \
                        sites found")

        self.add_option('no-spliced-alignment', bool, default=None,
                        optional=True, description="disable spliced alignment")

        # note to self just allow R F U and if paired in extend accordingly
        # otherwise mixed single paired will not work
        # Truseq is RF (R)
        self.add_option('rna-strandness', str, choices=["R", "F", "U"],
                        default=None, optional=False,
                        desription="Specify strand-specific information \
                        (unstranded); paired and are extended F->FR, R->RF")

        self.add_option('tmo', bool, default=None, optional=True,
                        description="Reports only those alignments within \
                        known transcriptome")

        self.add_option('dta', bool, default=None, optional=True,
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

        self.add_option('np', str, default=None, optional=True,
                        description="penalty for non-A/C/G/Ts in read/ref (1)")

        self.add_option('rdg', str, default=None, optional=True,
                        description="read gap open, extend penalties (5,3)")

        # TODO: rdf not in param list of hisat2, but rfg
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

        self.add_option('fr', bool, default=None, optional=False,
                        description="-1, -2 mates align fw/rev, rev/fw, \
                        fw/fw (--fr)")
        self.add_option('rf', bool, default=None, optional=False,
                        description="-1, -2 mates align fw/rev, rev/fw, \
                        fw/fw (--fr)")
        self.add_option('ff', bool, default=None, optional=False,
                        description="-1, -2 mates align fw/rev, rev/fw, \
                        fw/fw (--fr)")

        self.add_option('no-mixed', bool, default=None, optional=True,
                        description="suppress unpaired alignments for paired \
                        reads")

        self.add_option('no-discordant', bool, default=None, optional=True,
                        description="suppress discordant alignments for \
                        paired reads")

        # no-dovetail, no-contain, no-overlap?

        ###########
        # Output: #
        ###########

        # t, un, al, un-conc, al-conc?

        self.add_option('un-gz', bool, default=None, optional=True,
                        description="write unpaired reads that didn't align \
                        to <path>, gzip compress output")

        self.add_option('quiet', bool, default=None, optional=True,
                        description="print nothing to stderr except serious \
                        errors")

        # met-file, met-stderr, met?

        self.add_option('no-head', bool, default=None, optional=True,
                        description="supppress header lines, i.e. lines \
                        starting with @")

        self.add_option('no-sq', bool, default=None, optional=True,
                        description="supppress @SQ header lines")

        self.add_option('rg-id', str, default=None, optional=True,
                        description="puts sample name in rg")

        self.add_option('rg', str, default=None, optional=True,
                        description="add <text> ('lab:value') to @RG line of \
                        SAM header. (Note: @RG line only printed when --rg-id \
                        is set.)")

        self.add_option('omit-sec-seq', bool, default=None, optional=True,
                        description="put '*' in SEQ and QUAL fields for \
                        secondary alignments")

        ################
        # Performance: #
        ################

        ##########
        # Other: #
        ##########

        self.add_option('qc-filter', bool, default=None, optional=True,
                        description="filter out reads that are bad according \
                        to QSEQ filter")

        # seed?

        self.add_option('non-deterministic', bool, default=None, optional=True,
                        description="seed rand. gen. arbitrarily instead of \
                        using read attributes")

    def runs(self, run_ids_connections_files):
        flags = ["q", "qseq", "f", "c", "ignore-quals", "nofw", "dta",
                 "norc", "no-mixed",  "no-discordant", "quiet", "qc-filter",
                 "non-deterministic", "no-temp-splicesite",
                 "no-spliced-alignment", "tmo", "no-head", "no-sq",
                 "omit-sec-seq"]

        strflags = ["n-ceil", "ma", "mp", "sp", "np", "rdg", "score-min", "k",
                    "rg", "pen-cansplice", "pen-noncansplice",
                    "pen-canintronlen", "pen-noncanintronlen", "min-intronlen",
                    "max-intronlen", "known-splicesite-infile",
                    "novel-splicesite-outfile", "novel-splicesite-infile"]

        self.set_cores(self.get_option('cores'))

        res = [self.get_option('fr'),
               self.get_option('rf'),
               self.get_option('ff')]

        if sum(res) > 1:
            message = "too many stranded flags fr, rf, ff: %s"
            raise StandardError(message % (res))

        # Check if option values are valid
        if not os.path.exists(self.get_option('index') + '.1.ht2'):
            logger.error("Could not find index file: %s.*" %
                         self.get_option('index'))
            sys.exit(1)

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read'][0]
                sr_input = run_ids_connections_files[run_id]['in/second_read'][0]

                is_paired_end = True
                input_paths = [fr_input]

                if sr_input is None:
                    is_paired_end = False
                else:
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

                        for flag in strflags:
                            if self.is_option_set_in_config(flag):
                                hisat2.extend(['--' + flag,
                                              self.get_option(flag)])

                        # why -2?
                        hisat2.extend(['-p', str(self.get_option('cores') - 2),
                                       '-x', self.get_option('index')])

                        if is_paired_end:
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

                        metrics = run.add_output_file(
                            'metrics',
                            '%s-hisat2-metrics.txt' % run_id,
                            input_paths)
                        hisat2.extend(['--met-file', metrics])

                        if self.is_option_set_in_config('un-gz'):
                            if self.get_option('un-gz') is True:
                                unaligned = run.add_output_file(
                                    'unaligned',
                                    '%s-hisat2-unaligned.fastq.gz' % run_id,
                                    input_paths)
                                hisat2.extend(['--un-gz', unaligned])

                        hisat2_pipe.add_command(hisat2, stderr_path=log_stderr)
                        res = run.add_output_file(
                            'alignments',
                            '%s-hisat2-results.sam.gz' % run_id,
                            input_paths)

                        # Compress hisat2 output
                        pigz = [self.get_tool('pigz'),
                                '--stdout']
                        hisat2_pipe.add_command(pigz, stdout_path=res)
