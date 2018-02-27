import sys
import yaml
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class ChimPipe(AbstractStep):
    '''
    TODO

    ChimPipe is a tool to discover gene fusions in human paired-end RNA-Seq data.

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
        self.require_tool('chimpipe')
        self.require_tool('mkdir')
        self.require_tool('tar')
        self.require_tool('rm')

#adding options
        self.add_option('genome_index', str, optional=False, 
                        description="Reference genome index in GEM-format"
			)

        self.add_option('annotation', str, optional=False,
                        description="Reference gene annotation in gtf-format"
                        )

        self.add_option('transcriptome_index', str, optional=False,
                        description="Annotated transcriptome index in GEM-format"
                        )

        self.add_option('transcriptome_keys', str, optional=False,
                        description="Transcriptome to genome conversion keys"
                        )
        
        self.add_option('sample_ID', str, optional=False,
                        description="Identifier to be used in output file names"
                        )

        self.add_option('cores', str, default='6')
 
        self.add_option('consensus_seq', str, optional = True, 
			description="""Sequence pair of consensus splice site bases.TODO  """
			)
        
	self.add_option('library_type', str, optional = True, 
			description="TODO"
			)
	self.add_option('similarity', str, optional = True, 
			description="TODO"
			)


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


		#create folder structure

                my_temp_dir = run.get_output_directory_du_jour_placeholder()

                my_input  = os.path.join(my_temp_dir, 'input')
                my_output = os.path.join(my_temp_dir, 'output')
		my_cp_temp = os.path.join(my_temp_dir, 'run_temp')


		#create logfiles                
                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-chimpipe-log_stderr.txt' % run_id,
                    input_paths)

                log_stdout = run.add_output_file(
                    'log_stdout',
                    '%s-chimpipe-log_stdout.txt' % run_id,
                    input_paths)

              

                
                # init 
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pseudo_init:
			#create folders
                        make_dirs = [ self.get_tool('mkdir'), '-p',
                                      my_input,
                                      my_output,
				      my_cp_temp
                                      ]
                        pseudo_init.add_command(make_dirs)
			

               #Assemble chimpipe command
                with run.new_exec_group() as exec_group:
                    chimpipe = [
                        self.get_tool('chimpipe'),
                        '--fastq_1', fr_input,
                        '--fastq_2', sr_input,
                        '-g', self.get_option('genome_index'), 
                        '-a', self.get_option('annotation'),
			'-t', self.get_option('transcriptome_index'),
			'-k', self.get_option('transcriptome_keys'),
			'--sample-id', self.get_option('sample_ID'),
			'--threads', self.get_option('cores'),
			'--tmp-dir', my_cp_temp,
			'--output-dir', my_output,
			'--no-cleanup',
                    ]
	

	      	    if self.is_option_set_in_config('consensus_seq'):		       
	      	        if self.get_option('consensus_seq'):
	      	    	    chimpipe.extend([
	      	    	    '-C', self.get_option('consensus_seq'),
	      	    	    '-c', self.get_option('consensus_seq')])
	      
	      	    if self.is_option_set_in_config('library_type'):
	      	        if self.get_option('library_type'):
	      	    	    chimpipe.extend([
	      	    	    '-l', self.get_option('library_type')])
	      
	      	    if self.is_option_set_in_config('similarity'):
	      	        if self.get_option('similarity'):
	      	    	    chimpipe.extend([
	      	    	    '--similarity-gene-pairs', self.get_option('similarity')])
	

                    exec_group.add_command(chimpipe,  
                                       stderr_path=log_stderr, 
                                       stdout_path=log_stdout)
                        
                        
                        

 		with run.new_exec_group() as exec_group:
		
		#pack outfolder into tar/zip
			out_archive = run.add_output_file(
					'tar_archive',
					'%s-chimpipe-out.tar.gz'% run_id, input_paths)


			tar_output = [self.get_tool('tar'),
				'-czf', out_archive, 
				my_output ]

 			exec_group.add_command(tar_output)

		with run.new_exec_group() as exec_group:
                #remove temp dir

                        rm_temp = [self.get_tool('rm'),
                                '-r', my_input, my_output, my_cp_temp]

                        exec_group.add_command(rm_temp)
