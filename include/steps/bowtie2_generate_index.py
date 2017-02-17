import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Bowtie2GenerateIndex(AbstractStep):
    '''
    bowtie2-build builds a Bowtie index from a set of DNA sequences.
    bowtie2-build outputs a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2,
    .4.bt2, .rev.1.bt2, and .rev.2.bt2. In the case of a large index these
    suffixes will have a bt2l termination. These files together constitute the
    index: they are all that is needed to align reads to that reference.
    The original sequence FASTA files are no longer used by Bowtie 2 once the
    index is built.

    http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer

    typical command line::

        bowtie2-build [options]* <reference_in> <bt2_index_base>
    '''

    def __init__(self, pipeline):
        super(Bowtie2GenerateIndex, self).__init__(pipeline)
        self.set_cores(6)

        self.add_connection('in/reference_sequence')
        self.add_connection('out/bowtie_index')

        self.require_tool('dd')
        self.require_tool('pigz')
        self.require_tool('bowtie2-build')

        self.add_option('index-basename', str, optional = False,
                        description="Base name used for the bowtie2 index.")
        self.add_option('large-index', bool, optional = True,
                        description="Force bowtie2-build to build a large index,"
                        " even if the reference is less than ~ 4 billion "
                        "nucleotides long.")
        self.add_option('noauto', bool, optional = True,
                        description="Disable the default behavior whereby "
                        "bowtie2-build automatically selects values for the "
                        "--bmax, --dcv and --packed parameters according to "
                        "available memory. Instead, user may specify values for "
                        "those parameters. If memory is exhausted during "
                        "indexing, an error message will be printed; it is up "
                        "to the user to try new parameters.")
        self.add_option('packed', bool, optional = True,
                        description="Use a packed (2-bits-per-nucleotide) "
                        "representation for DNA strings. This saves memory but "
                        "makes indexing 2-3 times slower. Default: off. This is "
                        "configured automatically by default; use -a/--noauto "
                        "to configure manually.")
        self.add_option('bmax', int, optional = True,
                        description="The maximum number of suffixes allowed in "
                        "a block. Allowing more suffixes per block makes "
                        "indexing faster, but increases peak memory usage. "
                        "Setting this option overrides any previous setting for "
                        "--bmax, or --bmaxdivn. Default (in terms of the "
                        "--bmaxdivn parameter) is --bmaxdivn 4. This is "
                        "configured automatically by default; use -a/--noauto "
                        "to configure manually.")
        self.add_option('bmaxdivn', int, optional = True,
                        description="The maximum number of suffixes allowed in "
                        "a block, expressed as a fraction of the length of the "
                        "reference. Setting this option overrides any previous "
                        "setting for --bmax, or --bmaxdivn. Default: --bmaxdivn "
                        "4. This is configured automatically by default; use "
                        "-a/--noauto to configure manually.")
        self.add_option('dcv', int, optional = True,
                        description="Use <int> as the period for the "
                        "difference-cover sample. A larger period yields less "
                        "memory overhead, but may make suffix sorting slower, "
                        "especially if repeats are present. Must be a power of "
                        "2 no greater than 4096. Default: 1024. This is "
                        "configured automatically by default; use -a/--noauto "
                        "to configure manually.")
        self.add_option('nodc', bool, optional = True,
                        description="Disable use of the difference-cover "
                        "sample. Suffix sorting becomes quadratic-time in the "
                        "worst case (where the worst case is an extremely "
                        "repetitive reference). Default: off.")
        self.add_option('offrate', int, optional = True,
                        description="To map alignments back to positions on the "
                        "reference sequences, it's necessary to annotate "
                        "('mark') some or all of the Burrows-Wheeler rows with "
                        "their corresponding location on the genome. "
                        "-o/--offrate governs how many rows get marked: the "
                        "indexer will mark every 2^<int> rows. Marking more "
                        "rows makes reference-position lookups faster, but "
                        "requires more memory to hold the annotations at "
                        "runtime. The default is 5 (every 32nd row is marked; "
                        "for human genome, annotations occupy about 340 "
                        "megabytes).")
        self.add_option('ftabchars', int, optional = True,
                        description="The ftab is the lookup table used to "
                        "calculate an initial Burrows-Wheeler range with "
                        "respect to the first <int> characters of the query. "
                        "A larger <int> yields a larger lookup table but faster "
                        "query times. The ftab has size 4^(<int>+1) bytes. The "
                        "default setting is 10 (ftab is 4MB).")
        self.add_option('seed', int, optional = True,
                        description="Use <int> as the seed for pseudo-random "
                        "number generator.")
        self.add_option('cutoff', int, optional = True,
                        description="Index only the first <int> bases of the "
                        "reference sequences (cumulative across sequences) and "
                        "ignore the rest.")

        # Options for dd
        self.add_option('dd-blocksize', str, optional = True, default = "256k")

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['large-index', 'noauto', 'packed', 'bmax', 'bmaxdivn', 'dcv',
                   'nodc', 'offrate', 'ftabchars', 'seed', 'cutoff']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():
            # Get the basename
            index_basename = "%s-%s" % (
                self.get_option('index-basename'), run_id)

            with self.declare_run(index_basename) as run:
                with run.new_exec_group() as exec_group:
                    input_paths = run_ids_connections_files[run_id]\
                                  ['in/reference_sequence']

                    temp_files = list()

                    for input_path in input_paths:
                        # Is input gzipped fasta?
                        is_fasta_gz = False
                        is_fasta = False
                        if len([_ for _ in ['fa.gz', 'fasta.gz', 'fna.gz',
                                            'mfa.gz']\
                                if input_path.endswith(_)]) > 0:
                            is_fasta_gz = True
                        elif len([_ for _ in ['fa', 'fasta', 'fna', 'mfa'] \
                                  if input_path.endswith(_)]) > 0:
                            is_fasta = True

                        # Temporary file is always in FASTA format
                        temp_file = run.add_temporary_file(
                            prefix = "fifo-%s" %
                            os.path.splitext(os.path.basename(input_path))[0],
                            suffix = ".fa")

                        # 2. Output files to fifo
                        if is_fasta_gz:
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in chunks
                                dd_in = [self.get_tool('dd'),
                                         'ibs=%s' % self.get_option('dd-blocksize'),
                                         'if=%s' % input_path]
                                # 2.2 command: Uncompress data
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--stdout']
                                # 2.3 Write file chunks to fifo
                                dd_out = [self.get_tool('dd'),
                                          'obs=%s' % self.get_option('dd-blocksize'),
                                          'of=%s' % temp_file]
                                temp_files.append(temp_file)

                                unzip_pipe.add_command(dd_in)
                                unzip_pipe.add_command(pigz)
                                unzip_pipe.add_command(dd_out)
                        elif is_fasta:
                            dd = [self.get_tool('dd'),
                                  'bs=%s' % self.get_option('dd-blocksize'),
                                  'if=%s' % input_path,
                                  'of=%s' % temp_file]
                            temp_files.append(input_path)
                            exec_group.add_command(dd)
                        else:
                            logger.error("File %s does not end with any "
                                         "expected suffix (fastq.gz or "
                                         "fastq). Please fix that issue." %
                                         input_path)
                            sys.exit(1)
                    bowtie_build = [self.get_tool('bowtie2-build')]
                    # Add options
                    bowtie_build.extend(option_list)
                    # Add list of reference sequences
                    bowtie_build.append(','.join(temp_files))
                    # Add basename
                    bowtie_build.append(os.path.join(
                        run.get_output_directory_du_jour_placeholder(),
                        index_basename) )

                    exec_group.add_command(bowtie_build)

                    # Announce output files !!! Here might be dragons !!!
                    # 1.) input seqeunce length >4 billion leads to changed
                    #     file extension '.ebwtl'
                    # 2.) options --noref and --justref reduce number of
                    #     output files (therefore I disabled them)
                    # These issues need to be fixed or thought about!
                    run.add_output_file(
                        'bowtie_index',
                        '%s.1.bt2' % index_basename,
                        input_paths)
                    run.add_output_file(
                        'bowtie_index',
                        '%s.2.bt2' % index_basename,
                        input_paths)
                    run.add_output_file(
                        'bowtie_index',
                        '%s.3.bt2' % index_basename,
                        input_paths)
                    run.add_output_file(
                        'bowtie_index',
                        '%s.4.bt2' % index_basename,
                        input_paths)
                    run.add_output_file(
                        'bowtie_index',
                        '%s.rev.1.bt2' % index_basename,
                        input_paths)
                    run.add_output_file(
                        'bowtie_index',
                        '%s.rev.2.bt2' % index_basename,
                        input_paths)
