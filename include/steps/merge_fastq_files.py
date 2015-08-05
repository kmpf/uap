from abstract_step import *


class MergeFastqFiles(AbstractStep):
    '''
    Merge all .fastq(.gz) files of a sample.
    '''
    
    def __init__(self, pipeline):
        super(MergeFastqFiles, self).__init__(pipeline)
        
        self.set_cores(1) # muss auch in den Decorator
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        
        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('pigz')

    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            with self.new_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]

                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)
                    else:
                        suffix = ['.fastq', '.fastq.gz']
                        #input_files_by_suffix = dict()
                        #for s in suffix:
                        #input_files_by_suffix[s] = [x for x in input_paths
                        #                             if x.endswith(s)]

                        temp_files = list()
                        for input_path in input_paths:
                            # Gzipped files are unpacked first
                            # !!! Might be worth a try to use fifos instead of
                            #     temp files!!!
                            if input_path.endswith('fastq.gz'):
                                unzip_pipe = run.new_exec_group().add_pipeline()
                                temp_file = run.add_temporary_file(
                                    '%s' % os.path.basename(input_path) )
                                logger.info("Temp file: %s" % temp_file)
                                temp_files.append(temp_file)
                                # 1. command: Read file in 4MB chunks
                                dd1 = [self.get_tool('dd'),
                                      'ibs=4M',
                                      'if=%s' % input_path]
                                # 2. command: Uncompress file to STDOUT
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--stdout']
                                # 3. Write file in 4MB chunks
                                dd2 = [self.get_tool('dd'),
                                       'obs=4M',
                                       'of=%s' % temp_file]
                                
                                unzip_pipe.add_command(dd1)
                                unzip_pipe.add_command(pigz)
                                unzip_pipe.add_command(dd2)


                        cat_exec_group = run.new_exec_group()
                        cat = ['cat']
                        cat.extend(temp_files)
                    
                        cat_command = cat_exec_group.add_command(
                            cat,
                            stdout_path = run.add_output_file(
                                "%s" % read,
                                "%s%s.fastq" %
                                (run_id, read_types[read]),
                                input_paths))
