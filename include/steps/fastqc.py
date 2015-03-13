import sys
from abstract_step import *
import pipeline
import re
import process_pool
import shutil
import yaml
from yaml import dump 

class Fastqc(AbstractStep):
    '''
    | The fastqc step  is a wrapper for the fastqc tool. 
    | It generates some quality metrics for fastq files.
    | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
    | For this specific instance only the zip archive is preserved
    '''
    
    def __init__(self, pipeline):
        super(Fastqc, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read_fastqc_report')
        self.add_connection('out/first_read_log_stderr')
        self.add_connection('out/second_read_fastqc_report')
        self.add_connection('out/second_read_log_stderr')
        
        self.require_tool('fastqc')
        self.add_option('contaminent-file', str, optional =True)
        
    def declare_runs(self):
        found_files = dict()
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        paired_end_info = dict()
        # fetch all incoming run IDs which produce reads...
        for read in read_types.keys():
            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/%s' % read):
                if not run_id in found_files:
                    found_files[run_id] = dict()

                if not read in found_files[run_id]:
                    found_files[run_id][read] = list()
                found_files[run_id][read].extend(input_paths)

        for run_id in found_files.keys():
            with self.declare_run(run_id) as run:
                for read in found_files[run_id].keys():
                    #weired python way to get 'file' of 'file.bla.txt'
                    input_paths = found_files[run_id][read]
                    if input_paths == [None]:
                        run.add_empty_output_connection("%s_fastqc_report" % read)
                        run.add_empty_output_connection("%s_log_stderr" % read)
                    else:
                        input_base = os.path.basename(
                            input_paths[0]).split('.', 1)[0]
                    
                        # fastqc does not allow individual naming of files but 
                        # appends _fastqc to input file 
                        run.add_private_info(
                            '%s_fastqc_default_name' % read,
                            ''.join([input_base, '_fastqc']))
                        run.add_output_file(
                            "%s_fastqc_report" % read, 
                            "%s%s-fastqc.zip" % (run_id, read_types[read]),
                            input_paths)
                        run.add_output_file(
                            "%s_log_stderr" % read, 
                            "%s%s-fastqc-log_stderr.txt" % 
                            (run_id, read_types[read]),
                            input_paths)

    def execute(self, run_id, run):
        for read in ['first_read', 'second_read']:
            fastqc_out_dir = None
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    # Fastqc only allows to write to a directory
                    fastqc_out_dir =  self.get_output_directory_du_jour()

                    out_path = run.get_single_output_file_for_annotation(
                        '%s_fastqc_report' % read)

                    in_path  = run.get_input_files_for_output_file(out_path)

                    # set up processes
                    fastqc = [self.get_tool('fastqc'), 
                              '--noextract', '-o', fastqc_out_dir]
                    fastqc.extend(in_path)

                # create the pipeline and run it
                log_stderr = run.get_single_output_file_for_annotation(
                    '%s_log_stderr' % read)
                pipeline.append(fastqc, stderr_path = log_stderr)
                
            fastqc_default_name = run.get_private_info('%s_fastqc_default_name' % read)
            fastqc_report_basename  = fastqc_default_name + '.zip'
            full_path_zipped_fastqc_report = os.path.join(fastqc_out_dir,  fastqc_report_basename)
        
            try:
                os.rename(full_path_zipped_fastqc_report, out_path)
            except OSError:
                raise StandardError("os.rename failed of %s to %s" % full_path_zipped_fastqc_report, out_path)
