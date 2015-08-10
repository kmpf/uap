import sys
from abstract_step import *
import pipeline
import re
import yaml

class Cutadapt(AbstractStep):
    '''
    The cutadapt step can be used to clip adapter sequences from RNASeq reads.
    
    Any adapter may contain ``((INDEX))`` which will be replaced with every
    sample's index. The resulting adapter is checked for sanity and an
    exception is thrown if the adapter looks non-legit.
    '''
    
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)
        
        self.set_cores(3)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        self.add_connection('out/log_first_read')
        self.add_connection('out/log_second_read')
        
        self.require_tool('cat')
        self.require_tool('cutadapt')
        self.require_tool('dd')
        self.require_tool('fix_qnames')
        self.require_tool('mkfifo')
        self.require_tool('pigz')

        # Options for cutadapt
        self.add_option('adapter-type', str, optional=True, default='-a')
        self.add_option('adapter-R1', str, optional = True)
        self.add_option('adapter-R2', str, optional = True)
        self.add_option('adapter-file', str, optional=True)
        self.add_option('use_reverse_complement', bool, default=False)
        self.add_option('minimal-length', int, default=10, optional=True)

        self.add_option('fix_qnames', bool, default = False)

    def runs(self, run_ids_connections_files):
        '''
        
        '''
        ## Make sure the adapter type is one of -a, -b or -g 
        if self.is_option_set_in_config('adapter-type'):
            if not self.get_option('adapter-type') in set(['-a','-b','-g']):
                raise StandardError("Option 'adapter-type' must be "
                                    "either '-a', '-b', or '-g'!")


        read_types = {'first_read': 'R1', 'second_read': 'R2'}
        paired_end_info = dict()
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                for read in read_types:                
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]
                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)
                        run.add_empty_output_connection("log_%s" % read)
                    else:
                        paired_end_info[run_id] = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
                        # make sure that adapter-R1/adapter-R2 or adapter-file are
                        # correctly set 
                        # this kind of mutual exclusive option checking is a bit 
                        # tedious, so we do it here.
                        if paired_end_info[run_id]:
                            if ( not self.is_option_set_in_config('adapter-R2') and 
                                 not self.is_option_set_in_config('adapter-file') ):
                                raise StandardError(
                                    "Option 'adapter-R2' or " +
                                    "'adapter-file' required because " +
                                    "sample %s is paired end!" % run_id)
                        elif ( self.is_option_set_in_config('adapter-R2') and
                               not self.is_option_set_in_config('adapter-file') ):
                            raise StandardError(
                                "Option 'adapter-R2' not allowed because " +
                                "sample %s is not paired end!" % run_id)
                            
                        if ( self.is_option_set_in_config('adapter-file') and
                             self.is_option_set_in_config('adapter-R1') ):
                            raise StandardError(
                                "Option 'adapter-R1' and 'adapter-file' " +
                                "are both set but are mutually exclusive!" )
    
                        if ( not self.is_option_set_in_config('adapter-file') and
                             not self.is_option_set_in_config('adapter-R1') ):
                            raise StandardError(
                                "Option 'adapter-R1' or 'adapter-file' " +
                                "required to call cutadapt for sample %s!" % 
                                run_id)

                        temp_fifos = list()
                        exec_group = run.new_exec_group()
                        for input_path in input_paths:
                            # 1. Create temporary fifo for every input file
                            temp_fifo = run.add_temporary_file(
                                "fifo-%s" % os.path.basename(input_path) )
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)
                            # 2. Output files to fifo
                            if input_path.endswith('fastq.gz'):
                                with exec_group.add_pipeline() as pigz_pipe:
                                    # 2.1 command: Read file in 4MB chunks
                                    dd_in = [self.get_tool('dd'),
                                           'ibs=4M',
                                           'if=%s' % input_path]
                                    # 2.2 command: Uncompress file to fifo
                                    pigz = [self.get_tool('pigz'),
                                            '--decompress',
                                            '--stdout']
                                    # 2.3 command: Write file in 4MB chunks to 
                                    #              fifo
                                    dd_out = [self.get_tool('dd'),
                                              'obs=4M',
                                              'of=%s' % temp_fifo]

                                    pigz_pipe.add_command(dd_in)
                                    pigz_pipe.add_command(dd_pigz)
                                    pigz_pipe.add_command(dd_out)

                            elif input_path.endswith('fastq'):
                                # 2.1 command: Read file in 4MB chunks and
                                #              write to fifo in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                         'bs=4M',
                                         'if=%s' % input_path,
                                         'of=%s' % temp_fifo]
                                exec_group.add_command(dd_in)
                            else:
                                raise StandardError("File %s does not end with "
                                                    "any expected suffix ("
                                                    "fastq.gz or fastq). Please "
                                                    "fix that issue.")

                        # 3. Read data from fifos
                        with exec_group.add_pipeline() as cutadapt_pipe:
                            # 3.1 command: Read from ALL fifos
                            cat = [self.get_tool('cat')]
                            cat.extend(temp_fifos)
                            cutadapt_pipe.add_command(cat)
                            # 3.2 command: Fix qnames if user wants us to
                            if self.get_option('fix_qnames') == True:
                                fix_qnames = [self.get_tool('fix_qnames')]
                                cutadapt_pipe.add_command(fix_qnames)

                            # Let's get the correct adapter sequences or
                            # adapter sequence fasta file 
                            adapter = None
                            # Do we have adapter sequences as input?
                            if self.is_option_set_in_config('adapter-%s' \
                                                            % read_types[read]):
                                # Get adapter sequence 
                                adapter = self.get_option(
                                    'adapter-%s' % read_types[read])
                            
                                # add index to adapter sequence if necessary
                                if '((INDEX))' in adapter:
                                    index = self.find_upstream_info_for_input_paths(
                                        input_paths,
                                        'index-%s' % read_types[read])
                                    adapter = adapter.replace('((INDEX))', index)
                            
                                # create reverse complement if necessary
                                if self.get_option('use_reverse_complement'):
                                    complements = string.maketrans('acgtACGT',
                                                                   'tgcaTGCA')
                                    adapter = adapter.translate(complements)[::-1]
                                
                                # make sure the adapter is looking good
                                if re.search('^[ACGT]+$', adapter) == None:
                                    raise StandardError("Unable to come up with "
                                                        "a legit-looking adapter:"
                                                        "%s" % adapter)
                            # Or do we have a adapter sequence fasta file?
                            elif self.is_option_set_in_config('adapter-file'):
                                adapter = "file:" + self.get_option('adapter-file')

                            # 3.3 command: Clip adapters
                            cutadapt = [self.get_tool('cutadapt'), 
                                        self.get_option('adapter-type'), 
                                        adapter, '-']
                            cutadapt_log_file = run.add_output_file(
                                    'log_%s' % read,
                                    '%s-cutadapt-%s-log.txt'
                                    % (run_id, read_types[read]),
                                    input_paths)
                            # 3.4 command: Compress output
                            pigz = [self.get_tool('pigz'),
                                     '--blocksize', '4096',
                                     '--processes', '1',
                                     '--stdout']
                            # 3.5 command: Write to output file in 4MB chunks
                            clipped_fastq_file = run.add_output_file(
                                "%s" % read,
                                "%s_%s.fastq.gz" %
                                (run_id, read_types[read]),
                                input_paths)

                            dd = [self.get_tool('dd'),
                                  'obs=4M',
                                  'of=%s' % clipped_fastq_file]


                            cutadapt_pipe.add_command(cutadapt,
                                                      stderr_path =\
                                                      cutadapt_log_file)

                            cutadapt_pipe.add_command(pigz)
                            cutadapt_pipe.add_command(dd)
