import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')

class Hisat2(AbstractStep):
    '''
    HISAT2 is a fast and sensitive alignment program for mapping next-generation
    sequencing reads (both DNA and RNA) to a population of human genomes (as well
    as to a single reference genome).

    https://ccb.jhu.edu/software/hisat2/index.shtml
    '''

    def __init__(self, pipeline):

        super(Hisat2, self).__init__(pipeline)

        self.set_cores(6) # consider to adjust this default number

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/un')
        self.add_connection('out/al')
#        self.add_connection('out/unconc')
#        self.add_connection('out/alconc')
        self.add_connection('out/summary')
        self.add_connection('out/met')
        self.add_connection('out/sam')
        self.add_connection('out/log')
        self.add_connection('out/err')

        self.require_tool('hisat2')

        # main settings
        self.add_option('ht2-idx', str, optional=False,
                        description="Index filename prefix (minus trailing .X.ht2).")
        # <m[1|2]> and <sam> are in and output files that are set by uap

        # Input:
        self.add_option('q', bool, optional=True, # -
                        description='query input files are FASTQ .fq/.fastq (default)')
        self.add_option('qseq', bool, optional=True,
                        description='query input files are in Illumina\'s qseq format')
        self.add_option('f', bool, optional=True, # -
                        description='query input files are (multi-)FASTA .fa/.mfa')
        self.add_option('r', bool, optional=True, # -
                        description='query input files are raw one-sequence-per-line')
        self.add_option('c', bool, optional=True, # - .. do we need this option?
                        description='<m1>, <m2>, <r> are sequences themselves, not files')
        self.add_option('skip', int, optional=True,
                        description='skip the first <int> reads/pairs in the input (none)')
        self.add_option('upto', int, optional=True,
                        description='stop after first <int> reads/pairs (no limit)')
        self.add_option('trim5', int, optional=True,
                        description='trim <int> bases from 5\'/left end of reads (0)')
        self.add_option('trim3', int, optional=True,
                        description='trim <int> bases from 3\'/right end of reads (0)')
        self.add_option('phred33', bool, optional=True,
                        description='qualities are Phred+33 (default)')
        self.add_option('phred64', bool, optional=True,
                        description='qualities are Phred+64')
        self.add_option('int-quals', bool, optional=True,
                        description='qualities encoded as space-delimited integers')

        # Alignment:
        self.add_option('n-ceil', str, optional=True,
                        description='func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)')
        self.add_option('ignore-quals', bool, optional=True,
                        description='treat all quality values as 30 on Phred scale (off)')
        self.add_option('nofw', bool, optional=True,
                        description='do not align forward (original) version of read (off)')
        self.add_option('norc', bool, optional=True,
                        description='do not align reverse-complement version of read (off)')

        # Spliced Alignment:
        self.add_option('pen-cansplice', int, optional=True,
                        description='penalty for a canonical splice site (0)')
        self.add_option('pen-noncansplice', int, optional=True,
                        description='penalty for a non-canonical splice site (12)')
        self.add_option('pen-canintronlen', str, optional=True,
                        description='penalty for long introns (G,-8,1) with canonical splice sites')
        self.add_option('pen-noncanintronlen', str, optional=True,
                        description='penalty for long introns (G,-8,1) with noncanonical splice sites')
        self.add_option('min-intronlen', int, optional=True,
                        description='minimum intron length (20)')
        self.add_option('max-intronlen', int, optional=True,
                        description='maximum intron length (500000)')
        self.add_option('known-splicesite-infile', str, optional=True,
                        description='provide a list of known splice sites <path>')
        self.add_option('novel-splicesite-outfile', str, optional=True,
                        description='report a list of splice sites <path>')
        self.add_option('novel-splicesite-infile', str, optional=True,
                        description='provide a list of novel splice sites <path>')
        self.add_option('no-temp-splicesite', bool, optional=True,
                        description='disable the use of splice sites found')
        self.add_option('no-spliced-alignment', bool, optional=True,
                        description='disable spliced alignment')
        self.add_option('rna-strandness', str, optional=True,
                        choices=["unstranded", "FR", "RF"],
                        description='specify strand-specific information (unstranded)')
        self.add_option('tmo', bool, optional=True,
                        description='reports only those alignments within known transcriptome')
        self.add_option('dta', bool, optional=True,
                        description='reports alignments tailored for transcript assemblers')
        self.add_option('dta-cufflinks', bool, optional=True,
                        description='reports alignments tailored specifically for cufflinks')
        self.add_option('avoid-pseudogene', bool, optional=True,
                        description='tries to avoid aligning reads to pseudogenes (experimental'
                        'option).')
        self.add_option('no-templatelen-adjustment', bool, optional=True,
                        description='disables template length adjustment for RNA-seq reads')

        # Scoring:
        self.add_option('mp', str, optional=True,
                        description='max and min penalties for mismatch; '
                        'lower qual = lower penalty <6,2>')
        self.add_option('sp', str, optional=True,
                        description='max and min penalties for soft-clipping; '
                        'lower qual = lower penalty <2,1>')
        self.add_option('no-softclip', bool, optional=True,
                        description='no soft-clipping')
        self.add_option('np', int, optional=True,
                        description='penalty for non-A/C/G/Ts in read/ref (1)')
        self.add_option('rdg', str, optional=True,
                        description='read gap open, extend penalties (5,3)')
        self.add_option('rfg', str, optional=True,
                        description='--rfgreference gap open, extend penalties (5,3)')
        self.add_option('score-min', str, optional=True,
                        description='min acceptable alignment score w/r/t read length')

        # Reporting:
        self.add_option('k', int, optional=True, # -
                        description='(default: 5) report up to <int> alns per read')

        # Paired-end:
        self.add_option('minins', int, optional=True,
                        description='minimum fragment length (0), only valid with '
                        '--no-spliced-alignment')
        self.add_option('maxins', int, optional=True,
                        description='maximum fragment length (500), only valid with '
                        '--no-spliced-alignment')
        self.add_option('library_type', str, optional=True, # this must be parsed later!
                        choices=['fr','rf','ff'], default="fr",
                        description='-1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)')
        self.add_option('no-mixed', bool, optional=True,
                        description='suppress unpaired alignments for paired reads')
        self.add_option('no-discordant', bool, optional=True,
                        description='suppress discordant alignments for paired reads')

        # Output: # consider these as uap output connections!

        # hisat output option time is enabled by uap

        # hisat output options un, al, unc-conc, al-conc, summary-file, new-summary, met-file get
        # paths specified by uap

        # hisat output options quiet, met-stderr, met, no-head, no-sq, omit-sec-seq are ignored,
        # since uap wants to log everything possible

        # Performance:
        self.add_option('offrate', int, optional=True,
                        description='override offrate of index; must be >= index\'s offrate')
        self.add_option('threads', int, optional=True,
                        description='number of alignment threads to launch (1)')
        self.add_option('reorder', bool, optional=True,
                        description='force SAM output order to match order of input reads')
        self.add_option('mm', bool, optional=True,
                        description='use memory-mapped I/O for index; many \'hisat2\'s can share')

        # Other:
        self.add_option('qc-filter', bool, optional=True,
                        description='filter out reads that are bad according to QSEQ filter')
        self.add_option('seed', int, optional=True,
                        description='seed for random number generator (0)')
        self.add_option('non-deterministic', bool, optional=True,
                        description='seed rand. gen. arbitrarily instead of using read attributes')
        self.add_option('remove-chrname', bool, optional=True,
                        description='remove \'chr\' from reference names in alignment')
        self.add_option('add-chrname', bool, optional=True,
                        description='add \'chr\' to reference names in alignment ')



        # Options for 'dd' and 'pigz':
        self.add_option('dd-blocksize', str, optional = True, default = "2M")
        self.add_option('pigz-blocksize', str, optional = True, default = "2048")

    def runs(self, run_ids_connections_files):

        options = ['q', 'f', 'r', 'c', 'k']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()

        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))

        options = ['qseq', 'skip', 'upto', 'trim5', 'trim3', 'phred33', 'phred64',
                   'int-quals', 'n-ceil', 'ignore-quals', 'nofw', 'norc', 'pen-cansplice', 'pen-noncansplice',
                   'pen-canintronlen', 'pen-noncanintronlen', 'min-intronlen', 'max-intronlen',
                   'known-splicesite-infile', 'novel-splicesite-outfile', 'novel-splicesite-infile',
                   'no-temp-splicesite', 'no-spliced-alignment', 'rna-strandness', 'tmo', 'dta',
                   'dta-cufflinks', 'avoid-pseudogene', 'no-templatelen-adjustment', 'mp', 'sp', 'no-softclip',
                   'np', 'rdg', 'rfg', 'score-min', 'minins', 'maxins', 'no-mixed',
                   'no-discordant', 'offrate', 'threads', 'reorder', 'mm', 'qc-filter', 'seed',
                   'non-deterministic', 'remove-chrname', 'add-chrname']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        read_types = {'first_read': '_R1', 'second_read': '_R2'}

        # adjust cores to usr set number i.a.
        if 'threads' not in set_options:
            option_list.append('--threads')
            option_list.append(str(self.get_cores()))
        else:
            self.set_cores(self.get_option('threads'))

        # treat option library_type
        lt = self.get_option('library_type')
        option_list.append('--%s' % lt)

        # enable wallclock
        option_list.append('--time')


        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:

                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                input_paths = [ y for x in [fr_input, sr_input] \
                               for y in x if y != None ]

                # Do we have paired end data?
                # Do we need this information? It's implemented in library type, no?
                is_paired_end = True
                if sr_input == [None]:
                    is_paired_end = False

                with run.new_exec_group() as exec_group:

                    #  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
                    hisat2 = [self.get_tool('hisat2')]
                    hisat2.extend(option_list)

                    # outfiles that contain other information
                    un_outfile = run.add_output_file('un',
                                                     '%s-hisat2-unpaired_unaligned.fastq' % run_id,
                                                     input_paths)
                    al_outfile = run.add_output_file('al',
                                                     '%s-hisat2-unpaired_aligned.fastq' % run_id,
                                                     input_paths)
                    #un_conc_outfile = run.add_output_file('unconc',
                    #                                 '%s-hisat2-paired_unconcordantly.fastq.gz' % run_id,
                    #                                 input_paths)
                    #al_conc_outfile = run.add_output_file('alconc',
                    #                                 '%s-hisat2-paired_concordantly.fastq.gz' % run_id,
                    #                                 input_paths)
                    summary_outfile = run.add_output_file('summary',
                                                     '%s-hisat2-summary.txt' % run_id,
                                                     input_paths)
                    met_outfile = run.add_output_file('met',
                                                     '%s-hisat2-metrics.txt' % run_id,
                                                     input_paths)
                    hisat2.extend(['--un-gz', un_outfile,
                                   '--al-gz', al_outfile,
                                   '--summary-file', summary_outfile, '--new-summary',
                                   '--met-file', met_outfile])

                    # main input options
                    hisat2.extend(['-x', self.get_option('ht2-idx')])
                    if is_paired_end:
                        hisat2.extend(['-1', fr_input[0], '-2', sr_input[0]])
                    else:
                        hisat2.extend(['-U', fr_input[0]])

                    # main output options
                    sam_outfile = run.add_output_file('alignments',
                                                      '%s-hisat2-results.sam' % run_id,
                                                      input_paths)
                    hisat2.extend(['-S', sam_outfile])

                    # log files
                    ## wall clock and other info is provided here
                    log_outfile = run.add_output_file('log',
                                                      '%s-hisat2-stdout.txt' % run_id,
                                                      input_paths)
                    ## let's see what's provided here
                    err_outfile = run.add_output_file('err',
                                                      '%s-hisat2-stderr.txt' % run_id,
                                                      input_paths)

                    exec_group.add_command(hisat2,
                                           stdout_path = log_outfile,
                                           stderr_path = err_outfile)
