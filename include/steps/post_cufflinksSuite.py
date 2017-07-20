import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class Post_CufflinksSuite(AbstractStep):

    '''The cufflinks suite can be used to assembly new transcripts and
    merge those with known annotations. However, the output .gtf files
    need to be reformatted in several aspects afterwards. This step
    can be used to reformat and filter the cufflinksSuite .gtf file.
    '''

    def __init__(self, pipeline):
        super(Post_CufflinksSuite, self).__init__(pipeline)

        self.set_cores(6)

        # merged assembly 'merged.gft'
        self.add_connection('in/features') # combined.gtf
        # reformatted assembly
        self.add_connection('out/features') # filtered.gtf
        self.add_connection('out/log_stderr')

        self.require_tool('post_cufflinks_merge')
        self.require_tool('cat')

        self.add_option('remove_gencode', bool,
                        description='Hard removal of gtf line which match \'ENS\' in gene_name field',
                        default=False )
        self.add_option('remove_unstranded', bool,
                        description='Removes transcripts without strand specifity',
                        default=False)
        self.add_option('gene_name', str, optional=True,
                        description='String to match in gtf field gene_name for discarding',
                        default='ENS')
        self.add_option('remove_by_gene_name', bool,
                        description='Remove gtf if matches \'string\' in gene_name field',
                        default=False)
        # we may want to remove classcodes:
        # e,o,p,r,s
        self.add_option('class_list', str, optional=True,
                        description='Class codes to be removed; possible \'=,c,j,e,i,o,p,r,u,x,s,.\'',
                        default=None)
        self.add_option('filter_by_class', bool,
                        description='Remove gtf if any class is found in class_code field, requieres class_list',
                        default=False)
        # transport hyphenations to the final program call
        self.add_option('filter_by_class_and_gene_name', bool,
                        description='Combines remove-by-class and remove-by-gene-name',
                        default=False)

    def runs(self, run_ids_connections_files):
        
        # compile list of options
        options=['remove_gencode','remove_unstranded','gene_name','remove_by_gene_name',
                 'class_list','filter_by_class','filter_by_class_and_gene_name']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        run_id = "magic"

        with self.declare_run(run_id) as run:

            input_paths = run_ids_connections_files[run_id]['in/features']
            
            outfile = run.add_output_file('features', 
                                          '%s-filtered.gtf' % run_id, input_paths)
            logfile = run.add_output_file('log_stderr', 
                                          '%s-log_stderr.txt' % run_id, input_paths)
           	
           	
            # 1. create pipeline
            with run.new_exec_group() as as_exec_group:
           	
           	cat = [self.get_tool('cat'), input_paths]
           	post_cufflinks_merge = [self.get_tool('post_cufflinks_merge')]
           	post_cufflinks_merge.extend(option_list)
           	    
           	with as_exec_group.add_pipeline() as pipe:
           	    pipe.add_command(cat)
           	    pipe.add_command(post_cufflinks_merge,
           	                     stdout_path=outfile,
           	                     stderr_path=logfile)
           	        
