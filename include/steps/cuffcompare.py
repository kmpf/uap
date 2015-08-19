import sys
import re
import yaml

from ..abstract_step import *
from .. import pipeline
from .. import process_pool

class Cuffcompare(AbstractStep):
    
    def __init__(self, pipeline):
        super(Cuffcompare, self).__init__(pipeline)
        
        self.set_cores(1)

        self.add_connection('in/features')
        self.add_connection('out/features') # *.combined.gtf
        self.add_connection('out/loci')     # *.loci
        self.add_connection('out/stats')    # *.stats
        self.add_connection('out/tracking') # *.tracking
        self.add_connection('out/log_stderr')

        self.require_tool('cuffcompare')
        self.add_option('reference', str)


    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/features'):
            with self.declare_run(run_id) as run:

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one feature file.")

                run.add_output_file('features', '%s.combined.gtf' %run_id, input_paths)
                run.add_output_file('loci', '%s.loci' % run_id, input_paths)
                run.add_output_file('stats', '%s.stats' % run_id, input_paths)
                run.add_output_file('tracking', '%s.tracking' % run_id, input_paths)
                run.add_output_file('log_stderr', '%s-cuffcompare-log_stderr.txt' % run_id, input_paths)

                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                run.add_private_info('in-features', input_paths[0])

    def execute(self, run_id, run):

        cuffcompare_out_directory = self.get_temporary_path('cuffcompare-out')
        in_file = run.get_private_info('in-features')   

        with process_pool.ProcessPool(self) as pool:
            cuffcompare = [self.get_tool('cuffcompare'),
                           '-o', run_id,
                           '-r', self.get_option('reference'),
                           '-R', 
                           in_file]

            print(cuffcompare)
            try:
                os.mkdir(cuffcompare_out_directory)
            except OSError:
                pass
            os.chdir(cuffcompare_out_directory)

            pool.launch(cuffcompare, 
                        stderr_path = run.get_single_output_file_for_annotation('log_stderr')
                        )


        try:
            os.rename(os.path.join(cuffcompare_out_directory, '%s.combined.gtf' % run_id), 
                      run.get_single_output_file_for_annotation('features'))
        except OSError:
            raise StandardError('No file: %s' % os.path.join(cuffcompare_out_directory, 
                                                             '%s.combined.gtf' % run_id))

        try:
            os.rename(os.path.join(cuffcompare_out_directory, '%s.loci' % run_id), 
                      run.get_single_output_file_for_annotation('loci'))
        except OSError:
            raise StandardError('No file: %s' % os.path.join(cuffcompare_out_directory, 
                                                             '%s.loci' % run_id))

        try:
            os.rename(os.path.join(cuffcompare_out_directory, '%s.stats' % run_id), 
                      run.get_single_output_file_for_annotation('stats'))
        except OSError:
            raise StandardError('No file: %s' % os.path.join(cuffcompare_out_directory, 
                                                             '%s.stats' % run_id))

        try:
            os.rename(os.path.join(cuffcompare_out_directory, '%s.tracking' % run_id), 
                      run.get_single_output_file_for_annotation('tracking'))
        except OSError:
            raise StandardError('No file: %s' % os.path.join(cuffcompare_out_directory, 
                                                             '%s.tracking' % run_id))
