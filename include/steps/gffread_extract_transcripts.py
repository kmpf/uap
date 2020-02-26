from logging import getLogger
from abstract_step import AbstractStep
import os

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
        self.add_connection('in/anno')

        self.add_connection('out/fasta')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.require_tool('gffread')

        self.add_option('gtf', str, optional=True, default=None,
                        description="path to gtf file")

        self.add_option('output-fasta-name', str, optional=False, default=None,
                        description="name of the outputfile trancriptom myfasta.fa")

    def runs(self, run_ids_connections_files):
        run_id = self.get_option('output-fasta-name')
        # dependenency files (r for required)
        rfiles = []
        with self.declare_run(run_id) as run:
            cmd = [self.get_tool('gffread'), '-w']

            path = run.get_output_directory_du_jour_placeholder()
            cmd.append(path + '/' + self.get_option('output-fasta-name'))

            for __ , connection  in run_ids_connections_files.items():
                if 'in/fasta' in connection:
                    cmd.append('-g')
                    cmd.append(connection['in/fasta'][0])
                    rfiles.append(connection['in/fasta'][0])
                    continue


            if self.is_option_set_in_config('gtf'):
                cmd.append(os.path.abspath(self.get_option('gtf')))
            else:
                for __ , connection  in run_ids_connections_files.items():
                    if 'in/anno' in connection:
                        cmd.append(connection['in/anno'][0])
                        rfiles.append(connection['in/anno'][0])
                        continue



            stderr_file = "%s-gffread_extract_transcripts-log_stderr.txt" % (run_id)
            log_stderr = run.add_output_file("log_stderr",
                                             stderr_file, rfiles)
            stdout_file = "%s-gffread_extract_transcripts-log_stdout.txt" % (run_id)
            log_stdout = run.add_output_file("log_stdout",
                                             stdout_file, rfiles)

            run.add_output_file("fasta",
                                self.get_option('output-fasta-name'), rfiles)

            exec_group = run.new_exec_group()
            exec_group.add_command(cmd,
                                   stdout_path=log_stdout,
                                   stderr_path=log_stderr)
