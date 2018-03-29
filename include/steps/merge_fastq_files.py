import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class MergeFastqFiles(AbstractStep):
    '''
    This step concatenates all .fastq(.gz) files belonging to a certain sample.
    First and second read files are merged separately. The output files are
    gzipped.
    '''
    
    def __init__(self, pipeline):
        super(MergeFastqFiles, self).__init__(pipeline)
        
        self.set_cores(4) 
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        
        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('mkfifo')
        self.require_tool('pigz')

        #self.add_option('compress-output', bool, optional = True,
        #                default = True)

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "2M")
        self.add_option('pigz-blocksize', str, optional = True, default = "2048")

    def runs(self, run_ids_connections_files):

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]

                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)
                    else:
                        temp_fifos = list()
                        exec_group = run.new_exec_group()
                        for input_path in input_paths:
                            # Gzipped files are unpacked first
                            # !!! Might be worth a try to use fifos instead of
                            #     temp files!!!
                            # 1. Create temporary fifo
                            temp_fifo = run.add_temporary_file(
                                "fifo-%s" %
                                os.path.basename(input_path) )
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)

                            is_gzipped = True if os.path.splitext(input_path)[1]\
                                         in ['.gz', '.gzip'] else False

                            # 2. Output files to fifo
                            if is_gzipped:
                                with exec_group.add_pipeline() as unzip_pipe:
                                    # 2.1 command: Read file in 'dd-blocksize' chunks
                                    dd_in = [
                                        self.get_tool('dd'),
                                        'ibs=%s' %
                                        self.get_option('dd-blocksize'),
                                        'if=%s' % input_path
                                    ]
                                    # 2.2 command: Uncompress file to fifo
                                    pigz = [self.get_tool('pigz'),
                                            '--processes', str(self.get_cores()),
                                            '--decompress',
                                            '--blocksize', self.get_option('pigz-blocksize'),
                                            '--stdout']
                                    # 2.3 Write file in 'dd-blocksize' chunks to fifo
                                    dd_out = [
                                        self.get_tool('dd'),
                                        'obs=%s' %
                                        self.get_option('dd-blocksize'),
                                        'of=%s' % temp_fifo
                                    ]
                                
                                    unzip_pipe.add_command(dd_in)
                                    unzip_pipe.add_command(pigz)
                                    unzip_pipe.add_command(dd_out)
                            elif os.path.splitext(input_path)[1] in\
                                 ['.fastq', '.fq']:
                                # 2.1 command: Read file in 'dd-blocksize' chunks and
                                #              write to fifo in 'dd-blocksize' chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'bs=%s' %
                                    self.get_option('dd-blocksize'),
                                    'if=%s' % input_path,
                                    'of=%s' % temp_fifo
                                ]
                                exec_group.add_command(dd_in)
                            else:
                                logger.error("File %s does not end with any "
                                             "expected suffix (fastq.gz or "
                                             "fastq). Please fix that issue." %
                                             input_path)
                                sys.exit(1)
                        # 3. Read data from fifos
                        with exec_group.add_pipeline() as pigz_pipe:
                            # 3.1 command: Read from ALL fifos
                            cat = [self.get_tool('cat')]
                            cat.extend(temp_fifos)
                            pigz_pipe.add_command(cat)

                            # 3.2 Gzip output file
                            #if self.get_option('compress-output'):
                            pigz = [self.get_tool('pigz'),
                                    '--processes', str(self.get_cores()),
                                    '--blocksize', self.get_option('pigz-blocksize'),
                                    '--stdout']
                            pigz_pipe.add_command(pigz)

                            # 3.3 command: Write to output file in 'dd-blocksize' chunks
                            stdout_path = run.add_output_file(
                                "%s" % read,
                                "%s%s.fastq.gz" %
                                (run_id, read_types[read]),
                                input_paths)
                            dd = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize'),
                                'of=%s' % stdout_path
                            ]
                            pigz_pipe.add_command(dd)
