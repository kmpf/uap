from uaperrors import StepError
import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class FastxQualityStats(AbstractStep):
    '''
    fastx_quality_stats generates a text file containing quality information
    of the input fastq data.

    Documentation::

        http://hannonlab.cshl.edu/fastx_toolkit/

    Tested fastqc release: 0.0.13
    '''

    def __init__(self, pipeline):
        super(FastxQualityStats, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection(
            'in/second_read',
            optional = True)
        self.add_connection('out/first_read_quality_stats')
        self.add_connection(
            'out/second_read_quality_stats',
            optional = True)

        # Step was tested for cat (GNU coreutils) release 8.25
        self.require_tool('cat')
        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for FASTX Toolkit release 0.0.13
        self.require_tool('fastx_quality_stats')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

        # [Options for 'fastx_quality_stats':]
        self.add_option('new_output_format', bool, optional=True,
                        description="New output format (with more information "
                        "per nucleotide/cycle).")
        # no info in help for this option ?! -> excluded
        self.add_option('quality', int, default=33, optional=True)

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional=True, default="2M")
        self.add_option('pigz-blocksize', str, optional=True, default="2048")

    def runs(self, cc):
        # get a list of all read files we have to count
        sample_input_paths_dict = dict()
        reads_counts_files = dict()
        read_files = list()

        options = {'new_output_format': '-N', 'quality': '-Q'}
        option_list = list()
        for option in [o for o in options.keys()
                       if self.is_option_set_in_config(o)]:
            if isinstance(self.get_option(option), bool) and \
               self.get_option(option):
                option_list.append(options[option])
            else:
                option_list.append(options[option])
                option_list.append(str(self.get_option(option)))

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in cc.keys():
            cc.switch_run_id(run_id)
            with self.declare_run(run_id) as run:
                for read in read_types:
                    if not cc.exists_connection_for_run(f"in/{read}"):
                        continue
                    connection = 'in/%s' % read
                    input_paths = cc[run_id][connection]

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
                                os.path.basename(input_path))
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)

                            # 2. Output files to fifo
                            if input_path.endswith('fastq.gz'):
                                with exec_group.add_pipeline() as unzip_pipe:
                                    # 2.1 command: Read file in 'dd-blocksize'
                                    # chunks
                                    dd_in = [
                                        self.get_tool('dd'),
                                        'ibs=%s' %
                                        self.get_option('dd-blocksize'),
                                        'if=%s' % input_path
                                    ]
                                    # 2.2 command: Uncompress file to fifo
                                    pigz = [self.get_tool('pigz'),
                                            '--decompress',
                                            '--processes',
                                            str(self.get_cores()),
                                            '--blocksize',
                                            self.get_option('pigz-blocksize'),
                                            '--stdout']
                                    # 2.3 Write file in 'dd-blocksize' chunks
                                    # to fifo
                                    dd_out = [
                                        self.get_tool('dd'),
                                        'obs=%s' %
                                        self.get_option('dd-blocksize'),
                                        'of=%s' % temp_fifo
                                    ]

                                    unzip_pipe.add_command(dd_in)
                                    unzip_pipe.add_command(pigz)
                                    unzip_pipe.add_command(dd_out)
                            elif input_path.endswith('fastq'):
                                # 2.1 command: Read file in 'dd-blocksize' chunks and
                                # write to fifo in 'dd-blocksize' chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'bs=%s' % self.get_option('dd-blocksize'),
                                    'if=%s' % input_path,
                                    'of=%s' % temp_fifo
                                ]
                                exec_group.add_command(dd_in)
                            else:
                                raise StepError(
                                    self, "File %s does not end with any "
                                    "expected suffix (fastq.gz or "
                                    "fastq). Please fix that issue.")
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
                                                   stdout_path=fastx_qs_file)
