from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class KallistoIndex(AbstractStep):
    '''
    https://pachterlab.github.io/kallisto/manual
    Builds a kallisto index
    Usage: kallisto index [arguments] FASTA-files
    '''

    def __init__(self, pipeline):
        super(KallistoIndex, self).__init__(pipeline)

        self.set_cores(1)
        # input connections
        self.add_connection('in/fasta')


        # output connections
        self.add_connection('out/kallisto-index')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')


        # required tools
        self.require_tool('kallisto')

        # options
        self.add_option('index', str, optional=False, default=None,
                        description="Filename for the kallisto index")


        self.add_option('kmer-size', int, optional=True, default=31,
                        description="k-mer (odd) length (default: 31, max value: 31)")

        self.add_option('make-unique', bool, optional=True, default=None,
                        description="k-mer (odd) length (default: 31, max value: 31)")


    def runs(self, run_ids_connections_files):


        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                kallisto = [self.get_tool('kallisto'), 'index']

                path = run.get_output_directory_du_jour_placeholder()
                kallisto.extend(['-i', path + '/' + self.get_option('index')])

                kallisto.extend(['-k', str(self.get_option('kmer-size'))])

                if self.is_option_set_in_config('make-unique') and self.get_option('make-unique') :

                    kallisto.append('--make-unique')

                fasta_files = run_ids_connections_files[run_id]['in/fasta']
                if isinstance(fasta_files, list):
                    kallisto.extend(fasta_files)
                else:
                    kallisto.append(fasta_files)



                stderr_file = "%s-kallisto-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, fasta_files)
                stdout_file = "%s-kallisto-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, fasta_files)


                run.add_output_file("kallisto-index", 
                                    self.get_option('index'), fasta_files)


                kallisto_index = run.new_exec_group()
                kallisto_index.add_command(kallisto, 
                                        stdout_path=log_stdout,
                                        stderr_path=log_stderr)
