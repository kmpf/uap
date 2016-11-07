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
#        self.require_tool('cd')

        self.add_option('ref-gtf', str, optional=True,
                        description='A "reference" annotation GTF. The input assemblies are merged together with the reference GTF and included in the final output.')

    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():
            
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/features']
                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                    
                # check whether there's exactly one feature file
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one feature file.")

                in_file = input_paths[0]

                features_file = run.add_output_file('features',
                                                    '%s.combined.gtf' % run_id,
                                                    input_paths)
                loci_file     = run.add_output_file('loci',
                                                    '%s.loci' % run_id,
                                                    input_paths)
                stats_file    = run.add_output_file('stats',
                                                    '%s.stats' % run_id,
                                                    input_paths)
                tracking_file = run.add_output_file('tracking',
                                                    '%s.tracking' % run_id,
                                                    input_paths)
                log_err_file  = run.add_output_file('log_stderr',
                                                    '%s-cuffcompare-log_stderr.txt' % run_id,
                                                    input_paths)

                # replace cuffcompare with 'touch', check if the
                # outputfiles are still in the main dir
                cuffcompare = [self.get_tool('cuffcompare'), '-R', '-o', run_id,
                               '-r', str(self.get_option('ref-gtf')),
                               in_file]

                with run.new_exec_group() as cc_exec_group:
                    cc_exec_group.add_command(cuffcompare)
