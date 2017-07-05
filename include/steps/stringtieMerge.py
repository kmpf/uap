import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class StringtieMerge(AbstractStep):

    '''
    # stringtie --merge <gtf.list> > outputpat/outputname
    '''

    def __init__(self, pipeline):
        super(StringtieMerge, self).__init__(pipeline)

        self.set_cores(2)

        # all .gft assemblies from all samples that have been produced with stringtie
        self.add_connection('in/assembling')
        # merged assembly 'merged.gft'
        self.add_connection('out/assembling') # merged.gtf
        self.add_connection('out/assemblies') # input assemblies txt file
        self.add_connection('out/log_stderr')
        self.add_connection('out/run_log')

        self.require_tool('stringtie')
        self.require_tool('printf')
        self.require_tool('mkdir')
        self.require_tool('mv')


        self.add_option('G', str, optional=True,
                        description='reference annotation to include in the merging (GTF/GFF3)')
        self.add_option('m', int, optional=True,
                        description='minimum input transcript length to include in the merge (default: 50)')
        self.add_option('c', int, optional=True,
                        description='minimum input transcript coverage to include in the merge (default: 0)')
        self.add_option('F', float, optional=True,
                        description='minimum input transcript FPKM to include in the merge (default: 1.0)')
        self.add_option('T', float, optional=True,
                        description='minimum input transcript TPM to include in the merge (default: 1.0)')
        self.add_option('f', float, optional=True,
                        description='minimum isoform fraction (default: 0.01)')
        self.add_option('g', int, optional=True,
                        description='gap between transcripts to merge together (default: 250)')
        self.add_option('l', str, optional=True,
                        description='name prefix for output transcripts (default: MSTRG)')

        self.add_option('p', int, optional=True,
                        default=2, description='Number of cores')

        self.add_option('run_id', str, optional=True,
                        default="merge", description='uap specific sets runid')
        

    def runs(self, run_ids_connections_files):
       
        # reset cores to number of threads
        self.set_cores(self.get_option('p'))    


        # compile list of options
        options=['G', 'm', 'c', 'F', 'T', 'f', 'g', 'l']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))

        # get all paths to the stringtie assemblies from each sample
        stringtie_sample_gtf = []
        for run_id in run_ids_connections_files.keys():
            stringtie_sample_gtf.append(run_ids_connections_files[run_id]['in/assembling'][0])



        run_id = self.get_option('run_id')
        with self.declare_run(run_id) as run:
            
            # create the filename of the assemblies.txt file
            assemblies = [self.get_tool('printf'), '\n'.join(stringtie_sample_gtf)]
            # print assemblies
            
            assemblies_file = run.add_output_file('assemblies', 
                                                  '%s-stringtieMerge-assemblies.txt' % 
                                                  run_id, stringtie_sample_gtf)
           

            # 1. create assemblies file
            with run.new_exec_group() as exec_group:
                exec_group.add_command(assemblies, stdout_path = assemblies_file)
                with exec_group.add_pipeline() as stringtie_pipe:
                    res = run.add_output_file('assembling', 
                                              '%s-stringtieMerge-merged.gtf' % 
                                              run_id, stringtie_sample_gtf)

                    log_err_file = run.add_output_file('log_stderr', 
                                                       '%s-stringtieMerge-log_stderr.txt' % 
                                                       run_id, stringtie_sample_gtf)

                  

                    stringtieMerge = [self.get_tool('stringtie'), '--merge']
                    stringtieMerge.extend(option_list)
                    stringtieMerge.append(assemblies_file)                        
                    stringtie_pipe.add_command(stringtieMerge, 
                                               stderr_path = log_err_file, 
                                               stdout_path = res)

            
