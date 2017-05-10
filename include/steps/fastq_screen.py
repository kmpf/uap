import sys
import logging
from abstract_step import AbstractStep
import os

logger = logging.getLogger('uap_logger')

class FastqScreen(AbstractStep):
    '''
    fastq-tools random sample
     quick r1 work around
    '''
    
    def __init__(self, pipeline):
        super(FastqScreen, self).__init__(pipeline)
        cores = 10
        self.set_cores(cores)
    


        self.add_connection('in/first_read')
#        self.add_connection('in/second_read')
        self.add_connection('out/fqc_report')
        self.add_connection('out/log_stderr')
        self.add_connection('out/fqc_image')
        self.add_connection('out/first_read')


         
        # require_tool evtl. in abstract_step verstecken
        self.add_option('config', str , optional=False)
        self.add_option('cores', int, default=cores)
        self.add_option('nohits', bool, default=False)
        self.add_option('subset', int, default=False)

        self.require_tool('fastq_screen')
        self.require_tool('bowtie2')
        self.require_tool('mv')
        self.require_tool('rm')

    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        self.set_cores(self.get_option('cores'))

#        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        read_types = {'first_read': '_R1'}

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]
                    if input_paths == [None]:
                        run.add_empty_output_connection("%s_first_read" %
                                                        read)
                        run.add_empty_output_connection("%s_log_stderr" % read)
                    else:
                        for input_path in input_paths:
                            out_path       = run.add_output_file("fqc_report", 
                                                                 "%s%s_screen.txt" % (run_id,read_types[read]), [input_path])
                            out_path_image = run.add_output_file("fqc_image",  "%s%s_screen.png" % (run_id, read_types[read]), [input_path])
                            log_stderr = run.add_output_file("log_stderr", "%s%s-fastqscreen-log_stderr.txt" % (run_id, read_types[read]), [input_path])
                                                                                   

                            fastq_screen_exec_group  = run.new_exec_group()
                            fastq_screen = [self.get_tool('fastq_screen'), 
                                            '-conf', self.get_option('config')]

                            prefix  = "%s%s.fastq.gz" % (run_id,read_types[read])
                            fastq_screen.extend(['--optionalname', prefix])

                            if self.get_option('subset'): 
                                fastq_screen.extend(['--subset', str(self.get_option('subset'))])

                            if self.get_option('nohits'): 
                                nohits      = run.add_output_file("first_read", 
                                                                 "%s%s_no_hits.fastq.gz" % (run_id,read_types[read]), [input_path])
                                fastq_screen.extend(['--nohits', nohits])
                            else:
                                run.add_empty_output_connection("first_read")

                            fastq_screen.extend(['--outdir', run.get_output_directory_du_jour_placeholder(), input_path])
                            fastq_screen_exec_group.add_command(fastq_screen, stderr_path = log_stderr)
