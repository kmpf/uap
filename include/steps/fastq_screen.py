import sys
import logging
from abstract_step import AbstractStep
import os

logger = logging.getLogger('uap_logger')

class FastqScreen(AbstractStep):
    '''
    Fastq Screen (ver. 0.11.*)

    FastQ Screen allows you to screen a library of sequences in FastQ format
    against a set of sequence databases so you can see if the composition of
    the library matches with what you expect.

    '''

    def __init__(self, pipeline):
        super(FastqScreen, self).__init__(pipeline)
        cores = 10
        self.set_cores(cores)

        self.add_connection('in/first_read')
#        self.add_connection('in/second_read')
        self.add_connection('out/fqc_report')
        self.add_connection('out/fqc_image')
        self.add_connection('out/fqc_html')
        self.add_connection('out/first_read')
        self.add_connection('out/log_stdout')
        self.add_connection('out/log_stderr')

        self.add_option('config', str , optional=False,
                        description="Manually specify a location for the \
                        configuration.")

        self.add_option('cores', int, default=cores)

        self.add_option('nohits', bool, default=False,
                        description="Writes to a file the sequences that did \
                        not map to any of the specified genomes. This option \
                        is equivalent to specifying --tag --filter 0000 \
                        (number of zeros corresponds to the number of genomes \
                        screened).  By default the whole input file will be \
                        mapped, unless overridden by --subset.")

        self.add_option('subset', int, default=False,
                        description="Don't use the whole sequence file, but \
                        create a temporary dataset of this specified number \
                        of reads. The dataset created will be of approximately \
                        (within a factor of 2) of this size. If the real \
                        dataset is smaller than twice the specified size \
                        then the whole dataset will be used. Subsets will \
                        be taken evenly from throughout the whole original \
                        dataset. By Default FastQ Screen runs with this \
                        parameter set to 100000. To process an entire dataset \
                        however, adjust --subset to 0.")

        self.add_option('with_optionalname', bool, optional=False,
                        description="Indicates if the self implemented \
                        parameter --optionalname should be used. The \
                        parameter --optionalname cannot be modified \
                        and will only work with a patched version of \
                        fastq_screen. (see https://ribogit.izi.fraunhofer.de/oneButton/patched-uap-tools/fastq-screen)")

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
                            # prepare output files
                            file_pattern = "%s%s_screen.txt" % (run_id,read_types[read])
                            run.add_output_file("fqc_report", file_pattern, [input_path])

                            file_pattern = "%s%s_screen.png" % (run_id, read_types[read])
                            run.add_output_file("fqc_image", file_pattern, [input_path])


                            file_pattern = "%s%s_screen.html" % (run_id, read_types[read])
                            run.add_output_file("fqc_html", file_pattern, [input_path])


                            file_pattern = "%s%s-fastqscreen-log_stdout.txt" % (run_id, read_types[read])
                            log_stdout = run.add_output_file("log_stdout",
                                                             file_pattern,
                                                             [input_path])

                            file_pattern = "%s%s-fastqscreen-log_stderr.txt" % (run_id, read_types[read])
                            log_stderr = run.add_output_file("log_stderr",
                                                             file_pattern,
                                                             [input_path])

                            # build fastq_screen command
                            fastq_screen_exec_group  = run.new_exec_group()
                            fastq_screen = [self.get_tool('fastq_screen'),
                                            '-conf', self.get_option('config')]

                            # this parameter is only available in self patched fastq_screen version
                            if self.get_option('with_optionalname'):
                                prefix  = "%s%s.fastq.gz" % (run_id,read_types[read])
                                fastq_screen.extend(['--optionalname', prefix])

                            if self.get_option('subset'):
                                fastq_screen.extend(['--subset', str(self.get_option('subset'))])

                            if self.get_option('nohits'):
                                file_pattern = "%s%s_no_hits.fastq.gz" % (run_id,read_types[read])
                                nohits = run.add_output_file("first_read",
                                                             file_pattern,
                                                             [input_path])
                                fastq_screen.extend(['--nohits', nohits])
                            else:
                                run.add_empty_output_connection("first_read")

                            fastq_screen.extend(['--outdir', run.get_output_directory_du_jour_placeholder(), input_path])
                            fastq_screen_exec_group.add_command(fastq_screen, stdout_path = log_stdout, stderr_path = log_stderr)

