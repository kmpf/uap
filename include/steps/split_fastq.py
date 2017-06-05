from itertools import (takewhile,repeat)
import math

from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class SplitFastq(AbstractStep):
    '''

    '''

    def __init__(self, pipeline):
        super(SplitFastq, self).__init__(pipeline)

        # in connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        # out connections
        self.add_connection('out/out_first_read')
        # self.add_connection('out/out_second_read')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                    engine and threads ie")

        self.add_option('readcount', int, optional=False,
                        description="Number of reads per targetfile")

        # required tools
        self.require_tool('split_fastq')
        self.require_tool('wc')

    def get_line_count(self, filename):
        f = open(filename, 'rb')
        bufgen = takewhile(lambda x: x, (f.read(1024 * 1024) for _ in repeat(None)))
        return sum(buf.count(b'\n') for buf in bufgen)

    def runs(self, run_ids_connections_files):

        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():

            # TODO: create new run id

            with self.declare_run(run_id) as run:
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                readcount = self.get_option('readcount')

                # TODO: add r2 if exists
                input_fileset = [r1]

                # get lines to calculate output count
                line_count = self.get_line_count(r1)
                outfile_count = int(math.ceil(float(line_count/4)/readcount))

                # register output files
                for i in range(1, outfile_count + 1):
                    file_name = '%s_%s.fastq' % (run_id, i)
                    run.add_output_file('out_first_read', file_name, input_fileset)

                basename = run.get_output_directory_du_jour_placeholder() + '/' + run_id
                split_fastq = [self.get_tool('split_fastq'), '-i', r1, '-n', str(readcount), '-p', basename]
                sf_exec_group = run.new_exec_group()
                sf_exec_group.add_command(split_fastq)

                if run_ids_connections_files[run_id]['in/second_read'][0] == None:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]