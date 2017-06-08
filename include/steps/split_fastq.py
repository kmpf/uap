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
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')

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
        # TODO: if file is gzipped, first unzip!
        f = open(filename, 'rb')
        bufgen = takewhile(lambda x: x, (f.read(1024 * 1024) for _ in repeat(None)))
        return sum(buf.count(b'\n') for buf in bufgen)

    def runs(self, run_ids_connections_files):

        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():

            # TODO: set upper limit in yaml (file_count)
            index_list = range(1, 6)

            for index in index_list:
                new_run_id = '%s_%s' % (run_id, str(index))
                with self.declare_run(new_run_id) as run:
                    r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                    readcount = self.get_option('readcount')

                    input_fileset = [r1]
                    paired_end = False

                    # add r2 if exists
                    if run_ids_connections_files[run_id]['in/second_read'][0] != None:
                        r2 = run_ids_connections_files[run_id]['in/second_read'][0]
                        input_fileset.append(r2)
                        paired_end = True

                    # get lines to calculate output count
                    # TODO: take outfile count temporary from yaml param (file_count)
                    line_count = self.get_line_count(r1)
                    outfile_count = int(math.ceil(float(line_count/4)/readcount))

                    # register output files
                    for i in range(1, outfile_count + 1):
                        if i == index:
                            file_name = '%s_%s_r1.fastq' % (new_run_id, i)
                            run.add_output_file('first_read', file_name, input_fileset)

                            if paired_end:
                                file_name = '%s_%s_r2.fastq' % (new_run_id, i)
                                run.add_output_file('second_read', file_name, input_fileset)

                    basename = run.get_output_directory_du_jour_placeholder() + '/' + new_run_id
                    split_fastq_r1 = [
                        self.get_tool('split_fastq'),
                        '-i', r1,
                        '-n', str(readcount),
                        '-p', basename,
                        '-s', str(index),
                        '-m', 'r1'
                    ]
                    sf_exec_group = run.new_exec_group()
                    sf_exec_group.add_command(split_fastq_r1)

                    if paired_end:
                        split_fastq_r2 = [
                            self.get_tool('split_fastq'),
                            '-i', r2,
                            '-n', str(readcount),
                            '-p', basename,
                            '-s', str(index),
                            '-m', 'r2'
                        ]
                        sf_exec_group.add_command(split_fastq_r2)
