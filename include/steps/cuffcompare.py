import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class CuffCompare(AbstractStep):

    '''
    CuffCompare is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step compares a cufflinks assembly to 
    known annotation. For details about cuffcompare we refer to
    the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/

    '''

    def __init__(self, pipeline):
        super(CuffCompare, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/features')    # cuffmerge output
        self.add_connection('out/features')   # *.combined.gft
        self.add_connection('out/loci')       # *.loci
        self.add_connection('out/stats')      # *.stats
        self.add_connection('out/tracking')   # *.tracking
        self.add_connection('out/log_stderr')

        self.require_tool('cuffcompare')
        self.require_tool('mv')

        self.add_option('ref-gtf', str, optional=True,
                        description='A "reference" annotation GTF. The input assemblies are merged together with the reference GTF and included in the final output.')
        self.add_option('tmp_dir', str, optional=False,
                        description='A temporary directory for cuffcompare to write to.')

    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/features']
                
                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                    
                # check whether there's exactly one feature file
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one feature file.")

                # Check temporary directory
                if self.is_option_set_in_config('tmp_dir'):
                    if not os.path.isdir(self.get_option('tmp_dir')):
                        #dir not present
                        logger.error("Directory %s not found" % self.get_option('tmp_dir'))
                        sys.exit(1)
                    if not os.access(self.get_option('tmp_dir'), os.W_OK):
                        #not accessible
                        logger.error("Directory %s not accessible." % self.get_option('tmp_dir'))
                        sys.exit(1)

                features_file = run.add_output_file('features', '%s.combined.gtf' % run_id, input_paths)
                loci_file     = run.add_output_file('loci', '%s.loci' % run_id, input_paths)
                stats_file    = run.add_output_file('stats', '%s.stats' % run_id, input_paths)
                tracking_file = run.add_output_file('tracking', '%s.tracking' % run_id, input_paths)
                log_err_file  = run.add_output_file('log_stderr', '%s-cuffcompare-log_stderr.txt' % run_id, input_paths)

                cuffcompare = [self.get_tool('cuffcompare'), '-R', '-o', run_id,
                               '-r', str(self.get_option('ref-gtf')), 
                               input_paths[0]]
                # fixed (existing) temporary directory
                cuffcompare_out_path = str(self.get_option('tmp_dir'))
                
                # stderr.txt does not need to be moved, it is already in the temp folder of uap
                result_files = {
                    '%s.combined.gtf' % run_id : features_file,
                    '%s.loci' % run_id : loci_file,
                    '%s.stats' % run_id : stats_file,
                    '%s.tracking' % run_id : tracking_file
                }
            
#                print cuffcompare

                # 1. change to tmp dir and run cuffcompare
                with run.new_exec_group() as cc_exec_group:
                    # cuffcompare writes its output directly in the current
                    # working directory, which is known only at runtime
                    # therefore it needs a fixed (existing) temporary directory
                    # and we need to change to this directory before we run
                    # cuffcompare
                    os.chdir(cuffcompare_out_path) # change to tmp_dir
                    cc_exec_group.add_command(cuffcompare, stderr_path = log_err_file)

                # 2. mv output files from temp dir to final location            
                with run.new_exec_group() as mv_exec_group:
                    for orig, dest_path in result_files.iteritems():
                        orig_path = os.path.join(cuffcompare_out_path, orig)
                        mv = [self.get_tool('mv'), orig_path, dest_path]
                        mv_exec_group.add_command(mv)

