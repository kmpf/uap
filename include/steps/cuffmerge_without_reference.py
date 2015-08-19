import sys
import glob
import yaml
import hashlib

from ..abstract_step import *
from .. import misc
from .. import process_pool

class CuffmergeWithoutReference(AbstractStep):
    '''
    Cuffmerge is included in the cufflinks package.
    http://cufflinks.cbcb.umd.edu/manual.html
    '''
    
    def __init__(self, pipeline):
        super(CuffmergeWithoutReference, self).__init__(pipeline)
        self.set_cores(6)
        
        self.add_connection('in/features')
        self.add_connection('out/features')
        self.add_connection('out/assemblies_txt')
        self.add_connection('out/log_stderr')
        self.add_connection('out/run_log')

        self.require_tool('pigz')
        self.require_tool('cuffmerge')

        self.add_option('genome', str)
#        self.add_option('reference', str)
        self.add_option('run_id', str)
     


    def declare_runs(self):
        ### make sure files are available

        cufflinks_sample_gtf = []
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/features'):
            cufflinks_sample_gtf.append(input_paths[0])
            
        run_id = self.get_option('run_id')
        with self.declare_run(run_id) as run:
            run.add_output_file('features', '%s-cuffmerge-merged.gtf' % run_id, cufflinks_sample_gtf)
            run.add_output_file('run_log', '%s-cuffmerge-run.log' % run_id, cufflinks_sample_gtf)
            run.add_output_file('assemblies_txt', '%s-cuffmerge-assemblies.txt' % run_id, cufflinks_sample_gtf)
            run.add_output_file('log_stderr', '%s-cuffmerge-log_stderr.txt' % run_id, cufflinks_sample_gtf)
            run.add_private_info('cufflinks_sample_gtf', cufflinks_sample_gtf)
            

            
              
    def execute(self, run_id, run):
        
        assemblies  = self.get_temporary_path('assemblies.txt-', 'file')
        #cuffmerge_out_path = self.get_temporary_path('cuffmerge-out')
        cuffmerge_out_path = self.get_output_directory_du_jour()
        
        
        f  = open(assemblies, 'w')
        for i in run.get_private_info('cufflinks_sample_gtf'):
            f.write(i + '\n')
        f.close()
            
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cores = str(self.get_cores())
                cuffmerge = [self.get_tool('cuffmerge'),
                             '-o', cuffmerge_out_path,
                             '-s', self.get_option('genome'),
                             '-p', cores,
#                             '-g', self.get_option('reference'),
                             assemblies]
                log_stderr = run.get_single_output_file_for_annotation('log_stderr')
                pipeline.append(cuffmerge, stderr_path = log_stderr)                   


        os.rename(os.path.join(cuffmerge_out_path, 'merged.gtf'), run.get_single_output_file_for_annotation('features'))
        os.rename(assemblies, run.get_single_output_file_for_annotation('assemblies_txt'))
        os.rename(os.path.join(cuffmerge_out_path, 'logs/run.log'), run.get_single_output_file_for_annotation('run_log'))
        os.rmdir(os.path.join(cuffmerge_out_path, 'logs'))
        
