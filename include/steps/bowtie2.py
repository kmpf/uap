import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Bowtie2(AbstractStep):
    '''
    Bowtie2 is an ultrafast and memory-efficient tool for aligning sequencing
    reads to long reference sequences. It is particularly good at aligning reads
    of about 50 up to 100s or 1,000s of characters, and particularly good at
    aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the
    genome with an FM Index to keep its memory footprint small: for the human
    genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports
    gapped, local, and paired-end alignment modes.

    The input reads must come as .f[ast]q[.gz] files.

    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    typical command line::

        bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]

    This step wraps release: bowtie2
    '''

    def __init__(self, pipeline):
        super(Bowtie2, self).__init__(pipeline)

        self.set_cores(6)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')

        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')
        # Step was tested for bowtie2 release 2.2.9
        self.require_tool('bowtie2')

        # Options for bowtie2
        # 2B implemented: -q,-qseq,-f,-r,-c
        # -x
        self.add_option('index', str, optional=False,
                        description="Path to bowtie2 index (not containing file "
                        "suffixes).")
        # input
        self.add_option('skip', int, optional=True,
                        description="Skip the first <int> reads/pairs in the "
                        "input. Default: none")
        self.add_option('upto', int, optional=True,
                        description="Stop after the first <int> reads/pairs in "
                        "the input. Default: no limit.")
        self.add_option('trim5', int, optional=True,
                        description="Trim <int> bases from 5'/left end of reads "
                        "(default=0)")
        self.add_option('trim3', int, optional=True,
                        description="Trim <int> bases from 3'/right end of reads"
                        "(default=0)")
        self.add_option('phred33', bool, optional=True,
                        description="Qualities are Phred+33 (default)")
        self.add_option('phred64', bool, optional=True, default=False,
                        description="Qualities are Phred+64")
        self.add_option('int-quals', bool, optional=True,
                        description="Qualities encoded as space-delimited integers")
        # presets, for --end-to-end
        self.add_option('very-fast', bool, optional=True,
                        description="Preset, same as: -D 5 -R 1 -N 0 -L 22 -i S,0,2.50")
        self.add_option('fast', bool, optional=True,
                        description="Preset, same as: -D 10 -R 2 -N 0 -L 22 -i S,0,2.50")
        self.add_option('sensitive', bool, optional=True,
                        description="Preset, same as: -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 "
                        "(default)")
        self.add_option('very-sensitive', bool, optional=True,
                        description="Preset, same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50")
        # presets for --local
        self.add_option('very-fast-local', bool, optional=True,
                        description="Preset, same as: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00")
        self.add_option('fast-local', bool, optional=True,
                        description="Preset, same as: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75")
        self.add_option('sensitive-local', bool, optional=True,
                        description="Preset, same as: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 "
                        "(default)")
        self.add_option('very-sensitive-local', bool, optional=True,
                        description="Preset, same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50")

        # alignments
        self.add_option('N', int, optional=True,
                        description="Max # mismatches in seed alignment; can be 0 or 1 "
                        "(default=0)")
        self.add_option('L', int, optional=True,
                        description="length of seed substrings; must be >3, <32 "
                        "(default=22)")
        self.add_option('i', str, optional=True,
                        description="interval between seed substrings w/r/t read len "
                        "(default=\"S,1,1.15\")")
        self.add_option('n-ceil', str, optional=True,
                        description="func for max # non-A/C/G/Ts permitted in aln "
                        "(default=\"L,0,0.15\")")
        self.add_option('dpad', int, optional=True,
                        description="include <int> extra ref chars on sides of DP table "
                        "(default=15)")
        self.add_option('gbar', int, optional=True,
                        description="disallow gaps within <int> nucs of read extremes "
                        "(default=4)")
        self.add_option('ignore-quals', bool, optional=True,
                        description="treat all quality values as 30 on Phred scale")
        self.add_option('nofw', bool, optional=True,
                        description="do not align forward (original) version of read")
        self.add_option('norc', bool, optional=True,
                        description="do not align reverse-complement version of read")
        self.add_option('no-1mm-upfront', bool, optional=True,
                        description="do not allow 1 mismatch alignments before "
                        "attempting to scan for the optimal seeded alignments")
        self.add_option('end-to-end', bool, optional=True,
                        description="entire read must align; no clipping. Decide for"
                        "this or local option. (default)")
        self.add_option('local', bool, optional=True,
                        description="local alignment; ends might be soft clipped. "
                        "Decide for this or end-to-end option.")

        # scoring
        self.add_option('ma', int, optional=True,
                        description="match bonus (0 for end-to-end, 2 for local). "
                        "(default=0)")
        self.add_option('mp', int, optional=True,
                        description="max penalty for mismatch; lower qual = lower penalty "
                        "(default=6)")
        self.add_option('np', int, optional=True,
                        description="penalty for non-A/C/G/Ts in read/ref "
                        "(default=1)")
        self.add_option('rdg', str, optional=True,
                        description="read gap open, extend penalties "
                        "(default=\"5,3\")")
        self.add_option('rfg', str, optional=True,
                        description="reference gap open, extend penalties"
                        "(default=\"5,3\")")
        self.add_option('score-min', str, optional=True,
                        description="min acceptable alignment score w/r/t read length"
                        "(G,20,8 for local, L,-0.6,-0.6 for end-to-end)"
                        "(default=\"L,-0.6,-0.6\")")

        # reporting
        self.add_option('k', int, optional=True,
                        description="report up to <int> alns per read; MAPQ not meaningful")
        self.add_option('all', bool, optional=True,
                        description="report all alignments; very slow, MAPQ not meaningful")

        # effort
        self.add_option('D', int, optional=True,
                        description="give up extending after <int> failed extends in a row "
                        "(default=15)")
        self.add_option('R', int, optional=True,
                        description="for reads w/ repetitive seeds, try <int> sets of seeds "
                        "(default=2)")

        # Paired-end
        self.add_option('minins', int, optional=True,
                        description="minimum fragment length (default=0)")
        self.add_option('maxins', int, optional=True,
                        description="maximum fragment length (default=500)")
        self.add_option('fr', bool, optional=True,
                        description="-1, -2 mates align fw/rev (default)")
        self.add_option('rf', bool, optional=True,
                        description="-1, -2 mates align rev/fw")
        self.add_option('ff', bool, optional=True,
                        description="-1, -2 mates align fw/fw")
        self.add_option('no-mixed', bool, optional=True,
                        description="suppress unpaired alignments for paired reads")
        self.add_option('no-discordant', bool, optional=True,
                        description="suppress discordant alignments for paired reads")
        self.add_option('no-dovetail', bool, optional=True,
                        description="not concordant when mates extend past each other")
        self.add_option('no-contain', bool, optional=True,
                        description="not concordant when one mate alignment contains other")
        self.add_option('no-overlap', bool, optional=True,
                        description="not concordant when mates overlap at all")

        # output
        self.add_option('time', bool, optional=True,
                        description="print wall-clock time taken by search phases")
        # Options --un|al|un-conc|al-conc|un-gz|met-file
        # would produce output connections that are not default and known previously
        # This should be solved via temp. files
        # For now these options are ignored
        # which of these is written anyway??
        self.add_option('un', str, optional=True,
                        description="Path to write unpaired reads that didn't align")
        self.add_option('al', str, optional=True,
                        description="Path to write unpaired reads that aligned at least "
                        "once to")
        self.add_option('un-conc', str, optional=True,
                        description="Path to write pairs that didn't align concordantly")
        self.add_option('al-conc', str, optional=True,
                        description="write pairs that aligned concordantly at least once to")
        self.add_option('compress_add_output', bool, optional=True,
                        description="Ads -gz to the 4 options above to produce "
                        "compressed output.")
        self.add_option('quiet', bool, optional=True,
                        description="print nothing to stderr except serious errors")
        self.add_option('met-file', str, optional=True,
                        description="send metrics to file at")
        self.add_option('met-stderr', bool, optional=True,
                        description="send metrics to stderr")
        self.add_option('met', int, optional=True,
                        description="report internal counters & metrics every <int> secs "
                        "(default=1)")
        self.add_option('no-unal', bool, optional=True,
                        description="suppress SAM records for unaligned reads")
        self.add_option('no-head', bool, optional=True,
                        description="suppress header lines, i.e. lines starting with @")
        self.add_option('no-sq', bool, optional=True,
                        description="suppress @SQ header lines")
        self.add_option('rg-id', str, optional=True,
                        description="set read group id, reflected in @RG line and RG:Z: "
                        "opt field")
        self.add_option('rg', str, optional=True,
                        description="add <text> (lab:value) to @RG line of SAM header."
                        "Note: @RG line only printed when --rg-id is set.")
        self.add_option('omit-sec-seq', bool, optional=True,
                        description="put * in SEQ and QUAL fields for secondary alignments")

        # performance
        # reset self.set_cores to this one!
        self.add_option('threads', int, optional=True,
                        description="number of alignment threads to launch "
                        "(default=1)")
        self.add_option('reorder', bool, optional=True,
                        description="force SAM output order to match order of input reads")
        self.add_option('mm', bool, optional=True,
                        description="use memory-mapped I/O for index; many 'bowtie's "
                        "can share")

        # other
        self.add_option('qc-filter', bool, optional=True,
                        description="filter out reads that are bad according to QSEQ filter")
        self.add_option('seed', int, optional=True,
                        description="seed for random number generator (default=0)")
        self.add_option('non-deterministic', bool, optional=True,
                        description="seed rand. gen. arbitrarily instead of using read "
                        "attributes")


        # Options for dd
        self.add_option('dd-blocksize', str, optional = True, default = "2M")
        self.add_option('pigz-blocksize', str, optional = True, default = "2048")

    def runs(self, run_ids_connections_files):

        # Check if option values are valid
        if not os.path.exists(self.get_option('index') + '.1.bt2'):
            logger.error("Could not find index file: %s.*" %
                         self.get_option('index'))
            sys.exit(1)

        # compile all options set
        ## 1st all options that are given via --
        options1 = ['phred33', 'phred64',
                    'int-quals', 'very-fast', 'fast', 'sensitive', 'very-sensitive',
                    'very-fast-local', 'fast-local', 'sensitive-local',
                    'very-sensitive-local', 'ignore-quals', 'nofw', 'norc',
                    'no-1mm-upfront', 'end-to-end', 'local', 'all',
                    'fr', 'rf', 'ff', 'no-mixed', 'no-discordant', 'no-dovetail',
                    'no-overlap', 'time', 'quiet', 'met-stderr', 'no-unal', 'no-head',
                    'no-sq', 'omit-sec-seq', 'reorder', 'mm', 'qc-filter',
                    'non-deterministic',
                    # all non-bolean options:
                    'skip', 'upto', 'trim5', 'trim3', 'n-ceil',
                    'dpad', 'gbar', 'ma', 'mp', 'np', 'rdg', 'rfg', 'score-min',
                    'minins', 'maxins', 'un', 'al', 'un-conc', 'al-conc',
                    'met-file', 'met', 'rg-id', 'rg', 'threads', 'seed']


        ## 2nd all options that require a value, given with -
        options2 = ['N', 'L', 'i', 'k', 'D', 'R']

        # this will hold all bowtie options
        option_list = list()

        set_options1 = [option for option in options1 if \
                        self.is_option_set_in_config(option)]

        # collect set -- options
        for option in set_options1:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        # threads option can overwrite default # of cores for bowtie
        # and the cores variable
        if 'threads' not in set_options1:
            option_list.append('--threads')
            option_list.append(str(self.get_cores()))
        else:
            self.set_cores(self.get_option('threads'))

        # collect set - options
        set_options2 = [option for option in options2 if \
                        self.is_option_set_in_config(option)]

        for option in set_options2:
            option_list.append('-%s' % option)
            option_list.append(str(self.get_option(option)))



        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                input_paths = [ y for x in [fr_input, sr_input] \
                               for y in x if y !=None ]

                # Do we have paired end data and is it exactly one ?
                is_paired_end = True
                if sr_input == [None]:
                    is_paired_end = False

                # Tophat is run in this exec group
                with run.new_exec_group() as exec_group:
                    # Lists of fifos
                    fr_temp_fifos = list()
                    sr_temp_fifos = list()
                    # 1.

                    def prepare_input(input_path, exec_group, temp_fifos):
                        # Create temporary fifo
                        temp_fifo = run.add_temporary_file(
                            'in-fifo-%s' %
                            os.path.basename(input_path) )
                        mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                        exec_group.add_command(mkfifo)
                        temp_fifos.append(temp_fifo)
                        # Is input gzipped fasta?
                        is_fastq_gz = False
                        if len([_ for _ in ['fq.gz', 'fastq.gz']\
                                if input_path.endswith(_)]) == 1:
                            is_fastq_gz = True
                        # If yes we need to decompress it
                        if is_fastq_gz:
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in 'dd-blocksize' chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'bs=%s' % self.get_option('dd-blocksize'),
                                    'if=%s' % input_path
                                ]
                                unzip_pipe.add_command(dd_in)
                                # 2.2 command: Uncompress data
                                pigz = [self.get_tool('pigz'),
                                        '--processes', str(self.get_cores()),
                                        '--decompress',
                                        '--blocksize', self.get_option('pigz-blocksize'),
                                        '--stdout']
                                unzip_pipe.add_command(pigz)
                                # 2.3 Write file in 'dd-blocksize' chunks to fifo
                                dd_out = [
                                    self.get_tool('dd'),
                                    'obs=%s' % self.get_option('dd-blocksize'),
                                    'of=%s' % temp_fifo
                                ]
                                unzip_pipe.add_command(dd_out)
                        else:
                            dd = [
                                self.get_tool('dd'),
                                'bs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_path,
                                'of=%s' % temp_fifo
                            ]
                            exec_group.add_command(dd)

                        return (exec_group, temp_fifos)

                    for input_path in fr_input:
                        exec_group, fr_temp_fifos = prepare_input(
                            input_path, exec_group, fr_temp_fifos)
                    # And if we handle paired end data
                    if is_paired_end:
                        for input_path in sr_input:
                            exec_group, sr_temp_fifos = prepare_input(
                                input_path, exec_group, sr_temp_fifos)

                    # 3. Map reads using bowtie2
                    with exec_group.add_pipeline() as bowtie2_pipe:
                        # Assemble bowtie2 command
                        bowtie2 = [self.get_tool('bowtie2')]

                        bowtie2.extend(option_list)
                        bowtie2.extend(['-x', self.get_option('index')])

                        if is_paired_end:
                            bowtie2.extend([
                                '-1', ','.join(fr_temp_fifos),
                                '-2', ','.join(sr_temp_fifos)])
                        else:
                            bowtie2.extend(['-U', ','.join(fr_temp_fifos)])

                        bowtie2_pipe.add_command(bowtie2)
                        # Compress bowtie2 output
                        pigz = [self.get_tool('pigz'),
                                '--processes', str(self.get_cores()),
                                '--blocksize', self.get_option('pigz-blocksize'),
                                '--stdout']
                        bowtie2_pipe.add_command(pigz)
                        # Write bowtie2 output to file
                        dd = [
                            self.get_tool('dd'),
                            'obs=%s' % self.get_option('dd-blocksize'),
                            'of=%s' %
                            run.add_output_file(
                                'alignments',
                                '%s-bowtie2-results.sam.gz' % run_id,
                                input_paths
                            )
                        ]
                        bowtie2_pipe.add_command(dd)
