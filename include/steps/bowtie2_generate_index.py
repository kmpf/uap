import os

from abstract_step import AbstractStep

class Bowtie2GenerateIndex(AbstractStep):
    '''
    bowtie-build builds a Bowtie index from a set of DNA sequences.
    bowtie-build outputs a set of 6 files with suffixes .1.ebwt, .2.ebwt,
    .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt.
    (If the total length of all the input sequences is greater than about
    4 billion, then the index files will end in ebwtl instead of ebwt.)
    These files together constitute the index: they are all that is needed to
    align reads to that reference.
    The original sequence files are no longer used by Bowtie once the index is
    built.
    
    http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer
    
    typical command line::

        bowtie-build [options]* <reference_in> <ebwt_base>
    '''
    
    def __init__(self, pipeline):
        super(Bowtie2GenerateIndex, self).__init__(pipeline)
        self.set_cores(6)

        self.add_connection('in/reference_sequence')
        self.add_connection('out/bowtie_index')

        self.require_tool('dd')
        self.require_tool('pigz')
        self.require_tool('bowtie2-build')

        self.add_option('index-basename', str, optional = False)
        self.add_option('large-index', bool, optional = True)
        self.add_option('noauto', bool, optional = True)
        self.add_option('packed', bool, optional = True)
        self.add_option('bmax', int, optional = True)
        self.add_option('bmaxdivn', int, optional = True)
        self.add_option('dcv', int, optional = True)
        self.add_option('nodc', bool, optional = True)
        self.add_option('offrate', int, optional = True)
        self.add_option('ftabchars', int, optional = True)
        self.add_option('seed', int, optional = True)
        self.add_option('cutoff', int, optional = True)

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['large-index', 'noauto', 'packed', 'bmax', 'bmaxdivn', 'dcv',
                   'nodc', 'offrate', 'ftabchars', 'seed', 'cutoff']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool) and \
               self.get_option(option):
                option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():
            # Get the basename
            index_basename = "%s-%s" % (self.get_option('index-basename'),
                                        run_id)

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
                                # 2.1 command: Read file in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                         'ibs=4M',
                                         'if=%s' % input_path]
                                # 2.2 command: Uncompress data
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--stdout']
                                # 2.3 Write file in 4MB chunks to fifo
                                dd_out = [self.get_tool('dd'),
                                          'obs=4M',
                                          'of=%s' % temp_file]
                                temp_files.append(temp_file)
                                
                                unzip_pipe.add_command(dd_in)
                                unzip_pipe.add_command(pigz)
                                unzip_pipe.add_command(dd_out)
                        elif is_fasta:
                            dd = [self.get_tool('dd'),
                                  'bs=4M',
                                  'if=%s' % input_path, 
                                  'of=%s' % temp_file]
                            temp_files.append(input_path)
                            exec_group.add_command(dd)
                        else:
                            raise StandardError("File %s does not end with "
                                                "any expected suffix ("
                                                "fastq.gz or fastq). Please "
                                                "fix that issue." %
                                                input_path)

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
