import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class MergeAssembly(AbstractStep):
    '''
    This step merges a single gtf file that has been produced by a previous step 
    (e.g. cuffmerge or cuffcompare) with a reference annotation. No lines
    are discarded. The two files are simply concatenated and subsequently 
    sorted by position.
    '''
    
    def __init__(self, pipeline):
        super(MergeAssembly, self).__init__(pipeline)
        
        self.set_cores(1) 
        
        self.add_connection('in/features')
        self.add_connection('out/features')
        self.add_connection('out/log_stderr')
        
        self.require_tool('cat')
        self.require_tool('sort')

        # [Options for the step:]
        self.add_option('reference', str, optional = False,
                        description = "The reference annotation file that should be merged.")
        # [Options for 'sort':]
        self.add_option('temp-sort-dir', str, optional = False, 
                        description = 'Intermediate sort files are stored into this directory. '
                        'Note that this directory needs to be present before running this step.')

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:

                input_paths = run_ids_connections_files[run_id]['in/features']

                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                    
                # check whether there's exactly one feature file
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one feature file.")

                # this is the file that is going to be merged with the file provided in option 'reference'
                in_file = input_paths[0]

                # check for the existence of temporary sort directory
                if not os.path.isdir(self.get_option('temp-sort-dir')):
                    #dir not present
                    logger.error("Directory %s not found" % self.get_option('temp-sort-dir'))
                    sys.exit(1)
                if not os.access(self.get_option('temp-sort-dir'), os.W_OK):
                    #not accessible
                    logger.error("Directory %s not accessible." % self.get_option('temp-sort-dir'))
                    sys.exit(1)

                # this is the prefix without directory for uap add output file
                prefix = '%s_cuffcompare' %  run_id

                features_file = run.add_output_file('features',
                                                    '%s.allTranscripts.gtf' % prefix,
                                                    input_paths)
                log_err_file  = run.add_output_file('log_stderr',
                                                    '%s-log_stderr.txt' % prefix,
                                                    input_paths)

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:

                        # 1. concatenate both files
                        cat = [self.get_tool('cat'), in_file, self.get_option('reference')]
                        pipe.add_command(cat)

                        # 2. sort concatenated file by position
                        sort = [self.get_tool('sort'), 
                                '-T', os.path.join(self.get_option('temp-sort-dir'), run_id),
                                '-k', '1,1',
                                '-k', '4g,4',
                                '-k', '5g,5',
                                '-'] # the input coming from the cat cmd in the pipe
                        pipe.add_command(sort,
                                        stdout_path = features_file,
                                        stderr_path = log_err_file)

