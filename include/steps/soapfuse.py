import sys
import yaml
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SOAPfuse(AbstractStep):
    '''
    # TODO caro
    bla bla bla

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
        self.add_connection('out/log_stderr')

        # sf 

        self.require_tool('soapfuse')

        self.add_option('es', int, optional=True, default=8,
                        description="The step you want to end at 1-9"
                        )


        self.add_option('c', int, optional=True, 
                        description="""soapfuse config:
                        in config the following variables will be overwritten"
                        path to index: DB_db_dir
                        path to soapfuse bin: PG_pg_dir
                        path to soapfuse source: PS_ps_dir 
                        suffix for fastq: PA_all_fq_postfix (ex: *fgastq.gz)
                        cores: PA_all_process_of_align_software
                        """       )

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


                # bowtie2 is run in this exec group
                with run.new_exec_group() as exec_group:


                    with exec_group.add_pipeline() as soapfuse_pipe:
                        # Assemble soapfuse command
                        soapfuse = [
                            self.get_tool('soapfuse'),
                            '-p', str(self.get_option('cores') - 2),
                        ]
                
                        
                        log_stderr = run.add_output_file(
                                'log_stderr',
                                '%s-soapfuse-log_stderr.txt' % run_id,
                                                               input_paths)
                        res = run.add_output_file(
                                'text',
                                '%s-text.foo' % run_id,
                                input_paths
                            )



                        soapfuse_pipe.add_command(soapfuse, stderr_path=log_stderr)

