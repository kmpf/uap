from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class GffreadExtractTranscripts(AbstractStep):
    '''
    extract transcripts from gtf 
    http://ccb.jhu.edu/software/stringtie/gff.shtml
    write a fasta file with spliced exons for each GFF transcript
    gffread -w transcripts.fa -g /path/to/genome.fa transcripts.gtf
    '''

    def __init__(self, pipeline):
        super(GffreadExtractTranscripts, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/fasta')

        self.add_connection('out/fasta')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.require_tool('gffread')

        self.add_option('gtf', str, optional=False, default=None,
                        description="path to gtf file")

        self.add_option('output-fasta-name', str, optional=False, default=None,
                        description="name of the outputfile trancriptom myfasta.fa")

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                cmd = [self.get_tool('gffread'), '-w']

                path = run.get_output_directory_du_jour_placeholder()
                cmd.append(path + '/' + self.get_option('output-fasta-name'))


                fasta_file = run_ids_connections_files[run_id]['in/fasta']
                cmd.append('-g')
                cmd.append(fasta_file[0])
                
                cmd.append(self.get_option('gtf'))

                stderr_file = "%s-gffread_extract_transcripts-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, fasta_file)
                stdout_file = "%s-gffread_extract_transcripts-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, fasta_file)

                run.add_output_file("fasta",
                                    self.get_option('output-fasta-name'), fasta_file)

                exec_group = run.new_exec_group()
                exec_group.add_command(cmd,
                                       stdout_path=log_stdout,
                                       stderr_path=log_stderr)
