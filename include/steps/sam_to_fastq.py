import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SamToFastq(AbstractStep):
    '''
    Tbla bla bla    http://www.htslib.org/doc/samtools.html
    '''

    def __init__(self, pipeline):
        super(SamToFastq, self).__init__(pipeline)
        
        self.set_cores(8)
        
        self.add_connection('in/alignments')
        self.add_connection('out/first_read')
        

        self.require_tool('samtools')
        self.require_tool('pigz')

        self.add_option('F', int, optional = True)
        self.add_option('addF', int, optional = True)
        self.add_option('f', int, optional = True)


    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    logger.error("Expected exactly one alignments file.")
                    StandardError()
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

               

                out = run.add_output_file(
                    "first_read",
                    "%s_%s-samto.fastq.gz" %  (run_id, 'R1'),
                    input_paths) 
 
                    
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
                        # 1.1 command: Uncompress file to no fucking fifo 
                        if is_gzipped:
                            pigz = [self.get_tool('pigz'),
                                    '--decompress',
                                    '--processes', '1',
                                    '--stdout']
                            pigz.extend( input_paths)
                            pipe.add_command(pigz)

                            # 2. command: Convert to fastq
                            samtools = [self.get_tool('samtools'), 'fastq']

                            if self.is_option_set_in_config('f'):
                                samtools.extend(['-f', str(self.get_option('f'))])

                            if self.is_option_set_in_config('F'):
                                samtools.extend(['-F', str(self.get_option('F'))])


                            if self.is_option_set_in_config('addF'):
                                samtools.extend(['-F', str(self.get_option('addF'))])
                                
                            samtools.append('-')    
                            pipe.add_command(samtools)

                            #3 save fastq file

                            pigzc= [self.get_tool('pigz'), '--processes', '2',
                                   '--fast', '-']


                            pipe.add_command(pigzc, stdout_path=out )
