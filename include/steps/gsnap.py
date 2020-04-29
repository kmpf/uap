from logging import getLogger
from abstract_step import AbstractStep
import os

logger = getLogger('uap_logger')


class Gsnap(AbstractStep):
    '''
    '''

    def __init__(self, pipeline):
        super(Gsnap, self).__init__(pipeline)

        # input connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        # output connections
        self.add_connection('out/alignments')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        # required tools
        self.require_tool('gsnap')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                                    engine and threads ie")

        self.add_option('D', str, optional=False, default=None,
                        description="Genome directory")

        self.add_option('d', str, optional=False, default=None,
                        description="Genome database")

        self.add_option('t', int, optional=True, default=1,
                        description="Number of worker threads")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_fileset = []
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                input_fileset.append(r1)

                r2 = None
                if 'in/second_read' in run_ids_connections_files[run_id]:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]
                    input_fileset.append(r2)

                gsnap = [self.get_tool('gsnap'), '--gunzip']
                gsnap.extend(['-D', os.path.abspath(self.get_option('D'))])
                gsnap.extend(['-d', os.path.abspath(self.get_option('d'))])
                if self.is_option_set_in_config('t'):
                    gsnap.extend(['-t', str(self.get_option('t'))])
                # Batch-Mode
                gsnap.extend(['-B', '5'])

                # kmer size to use in genome database
                gsnap.extend(['-k', '12'])

                # Look for novel splicing
                gsnap.extend(['-N', '1'])

                # Maximum number of mismatches allowed
                # (if not specified, then
                # defaults to the ultrafast level of ((readlength+index_interval-1)/kmer - 2))
                # (By default, the genome index interval is 3, but this can be changed
                # by providing a different value for -q to gmap_build when processing
                # the genome.)
                #star.extend(['-m', '8.0'])

                # Another format type, other than default.
                # Currently implemented: sam, m8 (BLAST tabular format)
                gsnap.extend(['-A', 'sam'])

                # output
                gsnap.extend(['-o', 'gsnap_out.sam'])
                gsnap.extend(input_fileset)

                run.add_output_file(
                    "alignments", "gsnap_out.sam", input_fileset)

                stderr_file = "%s-gsnap-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, input_fileset)
                stdout_file = "%s-gsnap-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, input_fileset)

                gsnap_eg = run.new_exec_group()
                gsnap_eg.add_command(gsnap, stdout_path=log_stdout,
                                     stderr_path=log_stderr)
