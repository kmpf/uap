import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class FastxQualityStats(AbstractStep):
    '''
    fastx_quality_stats generates a text file containing quality information
    of the input FASTQ data.

    Documentation::

        http://hannonlab.cshl.edu/fastx_toolkit/
    '''

    def __init__(self, pipeline):
        super(FastxQualityStats, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read_quality_stats')
        self.add_connection('out/second_read_quality_stats')

        self.add_option('new_output_format', bool, default= True, optional=True)
        self.add_option('quality', int, default=33, optional=True)

        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('mkfifo')
        self.require_tool('fastx_quality_stats')
        self.require_tool('pigz')

    def runs(self, run_ids_connections_files):
        # get a list of all read files we have to count
        sample_input_paths_dict = dict()
        reads_counts_files = dict()
        read_files = list()

        options = {'new_output_format': '-N', 'quality': '-Q'}
        option_list = list()
        for option in [o for o in options.keys() \
                       if self.is_option_set_in_config(o)]:
            if isinstance(self.get_option(option), bool) and \
               self.get_option(option):
                option_list.append(options[option])
            else:
                option_list.append(options[option])
                option_list.append(str(self.get_option(option)))

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]

                    # Check for empty connections
                    if input_paths == [None]:
                        run.add_empty_output_connection(
                            "%s_quality_stats" % read)
                    else:
                        temp_fifos = list()
                        exec_group = run.new_exec_group()
                        for input_path in input_paths:
                            temp_fifo = run.add_temporary_file(
                                "fifo-%s" %
                                os.path.basename(input_path) )
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)

                            # 2. Output files to fifo
                            if input_path.endswith('fastq.gz'):
                                with exec_group.add_pipeline() as unzip_pipe:
                                    # 2.1 command: Read file in 4MB chunks
                                    dd_in = [self.get_tool('dd'),
                                           'ibs=4M',
                                           'if=%s' % input_path]
                                    # 2.2 command: Uncompress file to fifo
                                    pigz = [self.get_tool('pigz'),
                                            '--decompress',
                                            '--stdout']
                                    # 2.3 Write file in 4MB chunks to fifo
                                    dd_out = [self.get_tool('dd'),
                                              'obs=4M',
                                              'of=%s' % temp_fifo]

                                    unzip_pipe.add_command(dd_in)
                                    unzip_pipe.add_command(pigz)
                                    unzip_pipe.add_command(dd_out)
                            elif input_path.endswith('fastq'):
                                # 2.1 command: Read file in 4MB chunks and
                                #              write to fifo in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                         'bs=4M',
                                         'if=%s' % input_path,
                                         'of=%s' % temp_fifo]
                                exec_group.add_command(dd_in)
                            else:
                                logger.error("File %s does not end with any "
                                             "expected suffix (fastq.gz or "
                                             "fastq). Please fix that issue.")
                                sys.exit(1)
                        # 3. Read data from fifos and check quality stats
                        with exec_group.add_pipeline() as fastx_pipe:
                            # 3.1 command: Read from ALL fifos
                            cat = [self.get_tool('cat')]
                            cat.extend(temp_fifos)
                            # 3.2 command: Compute quality statistics
                            fastx_qs_file = run.add_output_file(
                                "%s_quality_stats" % read,
                                "%s%s.fastq.quality.tsv" %
                                (run_id, read_types[read]),
                                input_paths)
                            fastx_qs = [self.get_tool('fastx_quality_stats')]
                            fastx_qs.extend(option_list)
                            fastx_pipe.add_command(cat)
                            fastx_pipe.add_command(fastx_qs,
                                                   stdout_path = fastx_qs_file)
