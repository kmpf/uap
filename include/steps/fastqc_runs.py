import sys
import logging
import re
import shutil
import yaml
from yaml import dump 
from functools import wraps

from ..abstract_step import *
from .. import command
from .. import pipeline
from .. import process_pool
from .. import run as run_module

logger = logging.getLogger('uap_logger')

class Fastqc(AbstractStep):
    '''
    The fastqc step  is a wrapper for the fastqc tool. 
    It generates some quality metrics for fastq files.
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
    For this specific instance only the zip archive is preserved
    '''
    
    def __init__(self, pipeline):
        super(Fastqc, self).__init__(pipeline)
        
        self.set_cores(1) # muss auch in den Decorator
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read_fastqc_report')
        self.add_connection('out/first_read_fastqc_report_webpage')
        self.add_connection('out/first_read_log_stderr')
        self.add_connection('out/second_read_fastqc_report')
        self.add_connection('out/second_read_fastqc_report_webpage')
        self.add_connection('out/second_read_log_stderr')
        
        # require_tool evtl. in abstract_step verstecken
        self.require_tool('fastqc')
        self.require_tool('mkdir')
        self.require_tool('mv')

    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in list(run_ids_connections_files.keys()):
            with self.declare_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]
                    if input_paths == [None]:
                        run.add_empty_output_connection("%s_fastqc_report" %
                                                        read)
                        run.add_empty_output_connection("%s_log_stderr" % read)
                    else:
                        for input_path in input_paths:
                            # Get base name of input file
                            input_base = os.path.basename(input_path)\
                                                .split('.', 1)[0]
                            # Create temporary output directory
                            temp_dir = run.add_temporary_directory(
                                "%s" % input_base )
                            mkdir_exec_group = run.new_exec_group()
                            mkdir = [self.get_tool('mkdir'), temp_dir]
                            mkdir_exec_group.add_command(mkdir)
                            # 1. Run fastqc for input file
                            fastqc_exec_group = run.new_exec_group()
                            fastqc = [self.get_tool('fastqc'),
                                      '--noextract', '-o',
                                      temp_dir]
                            fastqc.append(input_path)
                            fastqc_command = fastqc_exec_group.add_command(
                                fastqc,
                                stderr_path = run.add_output_file(
                                    "%s_log_stderr" % read,
                                    "%s%s-fastqc-log_stderr.txt" %
                                    (run_id, read_types[read]),
                                    [input_path]) )
                            # 2. Move fastqc results to final destination
                            mv_exec_group = run.new_exec_group()
                            mv1 = [self.get_tool('mv'),
                                  os.path.join( temp_dir,
                                                ''.join([input_base,
                                                         '_fastqc.zip'])),
                                  run.add_output_file(
                                      "%s_fastqc_report" % read,
                                      "%s%s-fastqc.zip" %
                                      (run_id, read_types[read]),
                                      [input_path])]
                            mv2 = [self.get_tool('mv'),
                                  os.path.join( temp_dir,
                                                ''.join([input_base,
                                                         '_fastqc.html'])),
                                  run.add_output_file(
                                      "%s_fastqc_report_webpage" % read,
                                      "%s%s-fastqc.html" %
                                      (run_id, read_types[read]),
                                      [input_path])]

                            mv_exec_group.add_command(mv1)
                            mv_exec_group.add_command(mv2)
