from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Pizzly(AbstractStep):

    '''
    Pizzly is a tool to discover gene fusions
    in human paired-end RNA-Seq data.

    Paper:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5209911/

    Manual including typical usage:
    http://chimpipe.readthedocs.io/en/latest/index.html
   '''

    def __init__(self, pipeline):
        super(Pizzly, self).__init__(pipeline)
        self.set_cores(6)

# adding input/output connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('in/fusion.txt')

        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        self.add_connection('out/tar_archive')

# adding required tools
        self.require_tool('pizzly')
        self.require_tool('mkdir')
        self.require_tool('tar')
        self.require_tool('rm')

# adding options
        self.add_option('k-size', str, optional=False,
                        description="Size of k-mer")

        self.add_option('gtf-file', str, optional=False,
                        description="Reference transcript annotation in .gtf")

        self.add_option('number_of_mm', str, optional=True,
                        description="Number of allowed mismatches")

        self.add_option('insert-size', str, optional=True,
                        description="Maximum insert size")

        self.add_option('fasta-file', str, optional=False,
                        description="Transcripts in .fa ")

        self.add_option('cores', str, default='6')

        self.add_option('output_prefix', str, optional=False,
                        description="Prefix for output")

        self.add_option('cache', str, optional=True,
                        description="Cached gtf file created in a pizzly run")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:

                # create folder structure

                # my_output = run.get_output_directory_du_jour_placeholder()
                my_input = run_ids_connections_files[run_id]['in/fusion.txt'][0]

                # create logfiles
                input_paths = []

                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-pizzly-log_stderr.txt' % run_id,
                    input_paths)

                log_stdout = run.add_output_file(
                    'log_stdout',
                    '%s-pizzly-log_stdout.txt' % run_id,
                    input_paths)

                # Assemble pizzly command
                with run.new_exec_group() as exec_group:
                    pizzly = [
                        self.get_tool('pizzly'),
                        '-k', self.get_option('k-size'),
                        '--gtf', self.get_option('gtf-file'),
                        '--fasta', self.get_option('fasta-file'),
                        '--output', self.get_option('output_prefix')]

                    if self.is_option_set_in_config('cache'):
                        if self.get_option('cache'):
                            pizzly.extend([
                                '-cache', self.get_option('cache')])

                    if self.is_option_set_in_config('insert-size'):
                        if self.get_option('insert-size'):
                            pizzly.extend([
                                '--insert-size',
                                self.get_option('insert-size')])

                    if self.is_option_set_in_config('number_of_mm'):
                        if self.get_option('number_of_mm'):
                            pizzly.extend([
                                '--align-score',
                                self.get_option('number_of_mm')])

                    pizzly.extend([my_input])
                    print my_input

                    exec_group.add_command(pizzly,
                                           stderr_path=log_stderr,
                                           stdout_path=log_stdout)

#                with run.new_exec_group() as exec_group:
#                    # pack outfolder into tar/zip
#                    out_archive = run.add_output_file(
#                        'tar_archive',
#                        '%s-pizzly-out.tar.gz' % run_id, input_paths)
#
#                    tar_output = [self.get_tool('tar'),
#                                  '-czf', out_archive,
#                                  my_output]
#
#                    exec_group.add_command(tar_output)
#
#                with run.new_exec_group() as exec_group:
#                    # remove temp dir
#                    rm_temp = [self.get_tool('rm'),
#                               '-r',
#                               my_input,
#                               my_output]
#
#                    exec_group.add_command(rm_temp)
