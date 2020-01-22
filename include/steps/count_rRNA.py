from uaperrors import UAPError
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
        self.add_connection('out/report_rRNA')
        

        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cut')
        self.require_tool('grep')
        self.require_tool('sort')
        self.require_tool('uniq')




    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise UAPError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

               

                out = run.add_output_file(
                    "report_rRNA",
                    "%s_%s-rRNA_count.txt" %  (run_id, 'R1'),
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
                            samtools = [self.get_tool('samtools'), 'view', '-S' ,'-']
                            pipe.add_command(samtools)

                            #3 save fastq file
                            cuta  = [self.get_tool('cut'), '-f', '2,3,4']
                            pipe.add_command(cuta)
                            cutb  = [self.get_tool('cut'), '-f', '1', '-d' , '|']
                            pipe.add_command(cutb)
                            grep = [self.get_tool('grep'), '-v', '*']
                            pipe.add_command(grep)

                            cutc  = [self.get_tool('cut'), '-f', '1', '-d' , '_']
                            pipe.add_command(cutc)
                            
                            sorta = [self.get_tool('sort')]

                            pipe.add_command(sorta)

                            uniq = [self.get_tool('uniq'), '-c']
                            pipe.add_command(uniq)
                            sortb = [self.get_tool('sort')] 
                            pipe.add_command(sortb, stdout_path=out )
                            



