from logging import getLogger
from abstract_step import AbstractStep
import os

logger = getLogger('uap_logger')


class Pizzly(AbstractStep):

    '''
    Pizzly is a tool to discover gene fusions
    in human paired-end RNA-Seq data.

    Paper:
    https://www.biorxiv.org/content/early/2017/07/20/166322

    Manual including typical usage:
    https://github.com/pmelsted/pizzly
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
        self.add_connection('out/fusion_json')
        self.add_connection('out/fusion_fasta')
        self.add_connection('out/unfiltered_json')
        self.add_connection('out/unfiltered_fasta')

# adding required tools
        self.require_tool('pizzly')

# adding options
        self.add_option('k-mer_size', str, optional=False,
                        description="Size of k-mer")

        self.add_option('gtf-file', str, optional=False,
                        description="Reference transcript annotation in .gtf")

        self.add_option('align-score', str, optional=True,
                        description="Number of allowed mismatches")

        self.add_option('insert-size', str, optional=True,
                        description="Maximum insert size")

        self.add_option('fasta-file', str, optional=False,
                        description="Transcripts in .fa ")

        self.add_option('cores', str, default='6')

        self.add_option('cache', str, optional=True,
                        description="Cached gtf file created in a pizzly run")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:

                # create folder structure

                my_input = run_ids_connections_files[run_id]['in/fusion.txt'][0]

                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-pizzly-log_stderr.txt' % run_id,
                    [my_input])

                log_stdout = run.add_output_file(
                    'log_stdout',
                    '%s-pizzly-log_stdout.txt' % run_id,
                    [my_input])

                info_json = run.add_output_file(
                    'fusion_json',
                    '%s.json' % run_id,
                    [my_input])

                info_json = run.add_output_file(
                    'fusion_fasta',
                    '%s.fusions.fasta' % run_id,
                    [my_input])

                info_json = run.add_output_file(
                    'unfiltered_json',
                    '%s.unfiltered.json' % run_id,
                    [my_input])

                info_json = run.add_output_file(
                    'unfiltered_fasta',
                    '%s.unfiltered.fusions.fasta' % run_id,
                    [my_input])

                # Assemble pizzly command
                with run.new_exec_group() as exec_group:
                    pizzly = [
                        self.get_tool('pizzly'),
                        '-k', self.get_option('k-mer_size'),
                        '--gtf', os.path.abspath(self.get_option('gtf-file')),
                        '--fasta', os.path.abspath(self.get_option('fasta-file')),
                        '--output', run_id]

                    if self.is_option_set_in_config('cache'):
                        pizzly.extend([
                            '--cache', self.get_option('cache')])

                    if self.is_option_set_in_config('insert-size'):
                        pizzly.extend([
                            '--insert-size',
                            self.get_option('insert-size')])

                    if self.is_option_set_in_config('align-score'):
                        pizzly.extend([
                            '--align-score',
                            self.get_option('align-score')])

                    pizzly.extend([my_input])


                    exec_group.add_command(pizzly,
                                           stderr_path=log_stderr,
                                           stdout_path=log_stdout)

