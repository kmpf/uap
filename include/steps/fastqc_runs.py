import sys
import re
import shutil
import yaml
from yaml import dump 
from functools import wraps

from abstract_step import *
import command
import pipeline
import process_pool
import run as run_module
        

class Fastqc(AbstractStep):
    '''
    | The fastqc step  is a wrapper for the fastqc tool. 
    | It generates some quality metrics for fastq files.
    | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
    | For this specific instance only the zip archive is preserved
    '''
    
    def __init__(self, pipeline):
        super(Fastqc, self).__init__(pipeline)
        
        self.set_cores(1) # muss auch in den Decorator
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read_fastqc_report')
        self.add_connection('out/first_read_log_stderr')
        self.add_connection('out/second_read_fastqc_report')
        self.add_connection('out/second_read_log_stderr')
        
        # require_tool evtl. in abstract_step verstecken
        self.require_tool('fastqc') 


    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            run = self.new_run(run_id)
            for read in read_types:
                connection = 'in/%s' % read
                input_paths = run_ids_connections_files[run_id][connection]
                if input_paths == [None]:
                    run.add_empty_output_connection("%s_fastqc_report" % read)
                    run.add_empty_output_connection("%s_log_stderr" % read)
                else:
                    fastqc_exec_group = run.new_exec_group()
                    fastqc = [self.get_tool('fastqc'),
                              '--noextract', '-o',
                              self.get_output_directory_du_jour_placeholder()]
                    fastqc.extend(input_paths)
                    
                    fastqc_command = fastqc_exec_group.new_command(
                        fastqc,
                        stderr_path = run.add_output_file(
                            "%s_log_stderr" % read, 
                            "%s%s-fastqc-log_stderr.txt" % 
                            (run_id, read_types[read]),
                            input_paths))
                    
                    mv_exec_group = run.new_exec_group()
                    input_base = os.path.basename(
                        input_paths[0]).split('.', 1)[0]
                    mv = [self.get_tool('mv'),
                          os.path.join(
                              self.get_output_directory_du_jour_placeholder(),
                              ''.join([input_base, '_fastqc.zip'])),
                          run.add_output_file("%s_fastqc_report" % read,
                            "%s%s-fastqc.zip" % (run_id, read_types[read]),
                            input_paths)]
                    mv_command = mv_exec_group.new_command(mv)
