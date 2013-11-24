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
        
        self.add_connection('in/reads')
        self.add_connection('out/fastqc_report')
        self.add_connection('out/log_stderr')
        
        self.require_tool('fastqc')
        self.add_option('contaminent-file', str, optional =True)
                      
        
    def declare_runs(self):
        # fetch all incoming run IDs which produce reads...
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')

            # decide which read type we'll handle based on whether this is
            # paired end or not
            read_types = ['-R1']
            if  is_paired_end:
                read_types.append('-R2')


            # put input files into R1/R2 bins (or one single R1 bin)
            input_path_bins = dict()
            for _ in read_types:
                input_path_bins[_] = list()

            for path in input_paths:
                which =  misc.assign_string(os.path.basename(path), read_types)
                input_path_bins[which].append(path) 
                   

            # now declare runs
            for which in read_types:
                with self.declare_run("%s%s" % (run_id, which)) as run:
                    my_path = input_path_bins[which]
                    #weired python way to get 'file' of 'file.bla.txt'
                    input_base = os.path.basename(my_path[0]).split('.', 1)[0]

                    #fastqc does not allow individual naming of files but appends _fastqc to input file 
                    run.add_private_info('fastqc_default_name' , ''.join([input_base, '_fastqc']))
                    run.add_output_file("fastqc_report", "%s%s-fastqc.zip" % (run_id, which), input_path_bins[which])
                    run.add_output_file("log_stderr", "%s%s-fastqc-log_stderr.txt" % (run_id, which), input_path_bins[which])

    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                # Fastqc only allows to write to a directory 
                fastqc_out_dir =  self.get_output_directory_du_jour()

                out_path = run.get_single_output_file_for_annotation('fastqc_report')
                in_path  = run.get_input_files_for_output_file(out_path)

                
                
                # set up processes                              
                fastqc = [self.get_tool('fastqc'), '--noextract', '-o', fastqc_out_dir]
                fastqc.extend(in_path)

                
                # create the pipeline and run it
                log_stderr = run.get_single_output_file_for_annotation('log_stderr')
                pipeline.append(fastqc, stderr_path = log_stderr)                   




        fastqc_default_name = run.get_private_info('fastqc_default_name')
        fastqc_report_basename  = fastqc_default_name + '.zip'

        full_path_zipped_fastqc_report = os.path.join(fastqc_out_dir,  fastqc_report_basename)


        
        try:
            os.rename(full_path_zipped_fastqc_report, out_path)
        except OSError:
            raise StandardError("os.rename failed of %s to %s" % full_path_zipped_fastqc_report, out_path) 


 
        # in case of:
        #unzipped_fastqc_report = os.path.join(fastqc_out_dir,  fastqc_default_name)
        #try:
        #    shutil.rmtree(unzipped_fastqc_report)
        #except OSError:
        #    raise StandardError('removing unzipped dir failes')



               

