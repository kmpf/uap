import sys
import yaml
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SOAPfuse(AbstractStep):
    '''
    # TODO caro

    SOAPfuse is a tool to discover 

    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    
    typical command line::

        bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]
    '''
    
    def __init__(self, pipeline):
        super(SOAPfuse, self).__init__(pipeline)
        self.set_cores(6)



        #config 
        #sampleliste

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/text')
        self.add_connection('out/sf_sample_list')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        # sf 

        self.require_tool('soapfuse')
        self.require_tool('cp')
        self.require_tool('mkdir')
        self.require_tool('ln')
        self.require_tool('echo')

        self.add_option('es', int, optional=True, default=8,
                        description="The step you want to end at 1-9"
                        )



        self.add_option('c', str, optional=False, 
                        description="""soapfuse config:
                        in config the following variables will be overwritten"
                        path to index: DB_db_dir
                        path to soapfuse bin: PG_pg_dir
                        path to soapfuse source: PS_ps_dir 
                        suffix for fastq: PA_all_fq_postfix (ex: *fastq.gz)
                        cores: PA_all_process_of_align_software
                        """       )


        self.add_option('path_to_index_dir', str, optional=False,
                        description="Sets 'DB_db_dir' in Saopfuse config"
                        )

    
        self.add_option('path_to_sf_bin_dir', str, optional=False,
                        description="Sets 'TODFODB_db_dir' in Saopfuse config"
                        )

        self.add_option('path_to_sf_source', str, optional=False,
                        description="Sets 'TODFODB_db_dir' in Saopfuse config"
                        )

        
        self.add_option('suffix_for_fq_file', str, optional=False,
                        description="Sets 'TODFODB_db_dir' in Saopfuse config"
                        )

        self.add_option('read_length', int, optional=False,
                        description="Sets 'TODFODB_db_dir' in Saopfuse config"
                        )



        self.add_option('cores', int, default=6)
        
    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))


        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read

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




                my_temp_dir = run.get_output_directory_du_jour_placeholder()

                my_input  = os.path.join(my_temp_dir, 'input')
                my_output = os.path.join(my_temp_dir, 'output')
                my_config = os.path.join(my_temp_dir, 'input',
                                         os.path.basename( self.get_option('c')))

                my_sample_dir = os.path.join(my_input, "A", "L")
                r1 = run_id + "_1." +  self.get_option('suffix_for_fq_file')
                r2 = run_id + "_2." +  self.get_option('suffix_for_fq_file')
                my_sample_r1 = os.path.join(my_sample_dir, r1)



                my_sample_r2 = os.path.join(my_sample_dir, r2)


                res = run.add_output_file(
                    'text',
                    '%s-text.foo' % run_id,
                    input_paths)


                sl = run.add_output_file(
                    'sf_sample_list',
                    '%s-sl.txt' % run_id,
                    input_paths)


                
                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-soapfuse-log_stderr.txt' % run_id,
                    input_paths)

                log_stdout = run.add_output_file(
                    'log_stdout',
                    '%s-soapfuse-log_stdout.txt' % run_id,
                    input_paths)

                

                
                # init 
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pseudo_init:

                        make_dirs = [ self.get_tool('mkdir'), '-p',
                                      my_input,
                                      my_sample_dir, 
                                      my_output
                                      ]
                        pseudo_init.add_command(make_dirs)

                        cp_cmd = [self.get_tool('cp'), 
                                  self.get_option('c'),
                                  my_config]
                        pseudo_init.add_command(cp_cmd)



                with run.new_exec_group() as exec_group:

                    ln_sample = [self.get_tool('ln'), '-s',
                                 fr_input, my_sample_r1]
                    exec_group.add_command(ln_sample)


                    ln_sample = [self.get_tool('ln'), '-s',
                                sr_input, my_sample_r2]
                    exec_group.add_command(ln_sample)



                

                with run.new_exec_group() as exec_group:
                    t = ['A', 'L', run_id, str(self.get_option('read_length'))]
                    sf_list = '\t'.join(t)
                    
                    
                    echo_sf_list = [self.get_tool('echo'), sf_list]
                    exec_group.add_command(echo_sf_list,
                                           stdout_path=sl)
                    

                    


                #replace variables in config
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as replace_vars:
                        cat_cmd = ['cat', my_config]
                        replace_vars.add_command(cat_cmd)
                    
                        text = ' '.join(['DB_db_dir', '=', self.get_option('path_to_index_dir').replace("/","\\/")
                        ])
                        sed_arg = 's/DB_db_dir.*/' + text + '/'
                        sed_cmd = ['sed', sed_arg]
                        replace_vars.add_command(sed_cmd)

                        text = ' '.join(['PG_pg_dir', '=', self.get_option('path_to_sf_bin_dir').replace("/","\\/")])
                        sed_arg = 's/PG_pg_dir.*/' + text + '/'
                        sed_cmd = ['sed', sed_arg]
                        replace_vars.add_command(sed_cmd)

                        text = ' '.join(['PS_ps_dir', '=', self.get_option('path_to_sf_source').replace("/","\\/")])
                        sed_arg = 's/PS_ps_dir.*/' + text + '/'
                        sed_cmd = ['sed', sed_arg]
                        replace_vars.add_command(sed_cmd)

                        text = ' '.join(['PA_all_fq_postfix', '=', self.get_option('suffix_for_fq_file')])
                        sed_arg = 's/PA_all_fq_postfix.*/' + text + '/'
                        sed_cmd = ['sed', sed_arg]
                        replace_vars.add_command(sed_cmd)

                        text = ' '.join(['PA_all_process_of_align_software', '=', str(self.get_option('cores'))])
                        sed_arg = 's/PA_all_process_of_align_software.*/' + text + '/'
                        sed_cmd = ['sed', sed_arg]
                        replace_vars.add_command(sed_cmd,
                                                 stdout_path=res)
                    

          
                #Assemble soapfuse command
                with run.new_exec_group() as exec_group:
                    soapfuse = [
                        self.get_tool('soapfuse'),
                        '-fd', my_input,
                        '-c', res,
                        '-l', sl, 
                        '-o', my_output,
                    ]

                    exec_group.add_command(soapfuse,  
                                       stderr_path=log_stderr, 
                                       stdout_path=log_stdout)
                        
                        
                        


                
                    



