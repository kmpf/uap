import sys
import yaml
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class ChimPipe(AbstractStep):
    '''
    TODO

    FusionCatcher is a tool to discover gene fusions in human paired-end RNA-Seq data.

    Paper: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-2-r12

    Manual including required folder structure and typical usage:

    https://sourceforge.net/p/soapfuse/wiki/Home/
    
   '''
    
    def __init__(self, pipeline):
        super(ChimPipe, self).__init__(pipeline)
        self.set_cores(6)

# adding input/output connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        self.add_connection('out/tar_archive')

# adding required tools
        self.require_tool('fusioncatcher')
        self.require_tool('mkdir')
        self.require_tool('tar')
        self.require_tool('rm')

#adding options
        self.add_option('index', str, optional=False, 
                        description="Path to index folder")

        self.add_option('cores', str, default='6')

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:

                # Get lst of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read'][0]
                sr_input = run_ids_connections_files[run_id]['in/second_read'][0]

                # Do we have paired end data and is it exactly one ?
                is_paired_end = True
                input_paths = [fr_input]

                if sr_input == None:
                    logger.error("Not paired end")
                    sys.exit(1)
                else:
                    input_paths.append(sr_input)

		# create folder structure

                my_temp_dir = run.get_output_directory_du_jour_placeholder()

                my_input  = os.path.join(my_temp_dir, 'input')
                my_output = os.path.join(my_temp_dir, 'output')

		# create logfiles                
                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-fusioncatcher-log_stderr.txt' % run_id,
                    input_paths)

                log_stdout = run.add_output_file(
                    'log_stdout',
                    '%s-fusioncatcher-log_stdout.txt' % run_id,
                    input_paths)
                
                # init 
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pseudo_init:
			#create folders
                        make_dirs = [self.get_tool('mkdir'), '-p',
                                      my_input,
                                      my_output]
                        pseudo_init.add_command(make_dirs)
			
               # Assemble chimpipe command
                with run.new_exec_group() as exec_group:
                    fusioncatcher = [
                        self.get_tool('fusioncatcher'),
                        '-d',self.get_option('index') ,
                        '-i', my_input,
                        '-o', my_output, 
			'--threads', self.get_option('cores'),
			'--keep-unmapped-read',
			'--skip-filter-adapter',
			'--extract-buffer-size=1000000000'
                    ]
	
                    exec_group.add_command(fusioncatcher,  
                                       stderr_path=log_stderr, 
                                       stdout_path=log_stdout)
                        
 		with run.new_exec_group() as exec_group:
		
		#pack outfolder into tar/zip
			out_archive = run.add_output_file(
					'tar_archive',
					'%s-fusioncatcher-out.tar.gz'% run_id, input_paths)

			tar_output = [self.get_tool('tar'),
				'-czf', out_archive, 
				my_output ]

 			exec_group.add_command(tar_output)

		with run.new_exec_group() as exec_group:
                #remove temp dir

                        rm_temp = [self.get_tool('rm'),
                                '-r', my_input, my_output]

                        exec_group.add_command(rm_temp)
