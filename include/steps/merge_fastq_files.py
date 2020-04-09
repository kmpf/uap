from uaperrors import StepError
import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class MergeFastqFiles(AbstractStep):
    '''
    This step concatenates all .fastq(.gz) files belonging to a certain sample.
    First and second read files are merged separately. The output files are
    gzipped.
    '''

    def __init__(self, pipeline):
        super(MergeFastqFiles, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection(
            'in/second_read',
            optional = True)
        self.add_connection('out/first_read')
        self.add_connection(
            'out/second_read',
            optional = True)

        # Step was tested for cat (GNU coreutils) release 8.25
        self.require_tool('cat')
        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional=True, default="2M")
        self.add_option('pigz-blocksize', str, optional=True, default="2048")

    def runs(self, cc):

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in cc.keys():
            cc.switch_run_id(run_id)
            with self.declare_run(run_id) as run:
                for read in read_types:
                    if not cc.connection_exists(connection = f'in/{read}'):
                        continue
                    connection = 'in/%s' % read
                    input_paths = cc[run_id][connection]

                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)
                    else:
                        temp_fifos = list()
                        exec_group = run.new_exec_group()
                        for input_path in input_paths:
                            # Gzipped files are unpacked first
                            # !!! Might be worth a try to use fifos instead of
                            #     temp files!!!
                            # 1. Create temporary fifo
                            temp_fifo = run.add_temporary_file(
                                "fifo-%s" %
                                os.path.basename(input_path))
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)

                            is_gzipped = True if os.path.splitext(input_path)[1]\
                                in ['.gz', '.gzip'] else False

                            # 2. Output files to fifo
                            if is_gzipped:
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
                                            '--processes',
                                            str(self.get_cores()),
                                            '--decompress',
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
                            elif os.path.splitext(input_path)[1] in\
                                    ['.fastq', '.fq']:
                                # 2.1 command: Read file in 'dd-blocksize' chunks and
                                # write to fifo in 'dd-blocksize' chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'bs=%s' %
                                    self.get_option('dd-blocksize'),
                                    'if=%s' % input_path,
                                    'of=%s' % temp_fifo
                                ]
                                exec_group.add_command(dd_in)
                            else:
                                raise StepError(
                                    self, "File %s does not end with any "
                                    "expected suffix (fastq.gz or "
                                    "fastq). Please fix that issue." %
                                    input_path)
                        # 3. Read data from fifos
                        with exec_group.add_pipeline() as pigz_pipe:
                            # 3.1 command: Read from ALL fifos
                            cat = [self.get_tool('cat')]
                            cat.extend(temp_fifos)
                            pigz_pipe.add_command(cat)

                            # 3.2 Gzip output file
                            # if self.get_option('compress-output'):
                            pigz = [self.get_tool('pigz'),
                                    '--processes',
                                    str(self.get_cores()),
                                    '--blocksize',
                                    self.get_option('pigz-blocksize'),
                                    '--stdout']
                            pigz_pipe.add_command(pigz)

                            # 3.3 command: Write to output file in
                            # 'dd-blocksize' chunks
                            stdout_path = run.add_output_file(
                                "%s" % read,
                                "%s%s.fastq.gz" %
                                (run_id, read_types[read]),
                                input_paths)
                            dd = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize'),
                                'of=%s' % stdout_path
                            ]
                            pigz_pipe.add_command(dd)
