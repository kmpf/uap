from uaperrors import StepError
import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class FixCutadapt(AbstractStep):
    '''
    This step takes FASTQ data and removes both reads of a paired-end read, if
    one of them has been completely removed by cutadapt (or any other software).
    '''

    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection(
            'in/second_read',
            optional = True)
        self.add_connection('out/first_read')
        self.add_connection(
            'out/second_read',
            optional = True)

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional=True, default="2M")
        self.add_option('pigz-blocksize', str, optional=True, default="2048")

        # Step was tested for cat (GNU coreutils) release 8.25
        self.require_tool('cat')
        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        self.require_tool('fix_cutadapt')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

    def runs(self, cc):

        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        for run_id in cc.keys():
            cc.switch_run_id(run_id)
            with self.declare_run(run_id) as run:
                temp_fifos = dict()
                exec_group = run.new_exec_group()
                for read in read_types:
                    if not cc.connection_exists(f"in/{read}"):
                        continue
                    connection = f"in/{read}"
                    input_paths = cc[run_id][connection]
                    temp_fifos[f"{read}_in"] = None
                    temp_fifos[f"{read}_out"] = None
                    if input_paths == [None]:
                        run.add_empty_output_connection(f"{read}")

                    elif len(input_paths) != 1:
                        raise StepError(
                            self, "Expected single input file. Found files "
                            "%s for run: %s" %
                            (input_paths, run_id))
                    else:
                        # 1. Create temporary fifos
                        # 1.1 Input fifo
                        temp_fifos[f"{read}_in"] = run.add_temporary_file(
                            "in-fifo-%s" %
                            os.path.basename(input_paths[0]))
                        mkfifo_in = [self.get_tool('mkfifo'),
                                     temp_fifos[f"{read}_in"]]
                        exec_group.add_command(mkfifo_in)
                        # 1.2 Output fifo
                        temp_fifos[f"{read}_out"] = run.add_temporary_file(
                            f"{read}-out-fifo")
                        mkfifo_out = [self.get_tool('mkfifo'),
                                      temp_fifos[f"{read}_out"]]
                        exec_group.add_command(mkfifo_out)

                        # 2. Output files to fifo
                        if input_paths[0].endswith('fastq.gz'):
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in 'dd-blocksize'
                                # chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'ibs=%s' % self.get_option('dd-blocksize'),
                                    'if=%s' % input_paths[0]
                                ]
                                # 2.2 command: Uncompress file to fifo
                                pigz = [self.get_tool('pigz'),
                                        '--processes',
                                        str(self.get_cores()),
                                        '--decompress',
                                        '--blocksize',
                                        self.get_option('pigz-blocksize'),
                                        '--stdout']
                                # 2.3 Write file in 'dd-blocksize' chunks to
                                # fifo
                                dd_out = [
                                    self.get_tool('dd'),
                                    'obs=%s' % self.get_option('dd-blocksize'),
                                    'of=%s' % temp_fifos["%s_in" % read]
                                ]

                                unzip_pipe.add_command(dd_in)
                                unzip_pipe.add_command(pigz)
                                unzip_pipe.add_command(dd_out)
                        elif input_paths[0].endswith('fastq'):
                            # 2.1 command: Read file in 'dd-blocksize' chunks and
                            # write to fifo in 'dd-blocksize' chunks
                            dd_in = [
                                self.get_tool('dd'),
                                'bs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_paths[0],
                                'of=%s' % temp_fifos["%s_in" % read]
                            ]
                            exec_group.add_command(dd_in)
                        else:
                            raise StepError(
                                self, "File %s does not end with any "
                                "expected suffix (fastq.gz or fastq). "
                                "Please fix that issue." %
                                input_paths[0])
                # 3. Start fix_cutadapt
                fix_cutadapt = [self.get_tool('fix_cutadapt'),
                                temp_fifos["first_read_in"],
                                temp_fifos["first_read_out"]]
                if cc.connection_exists("in/second_read"):
                    fix_cutadapt.extend([
                        '--R2-in', temp_fifos["second_read_in"],
                        '--R2-out', temp_fifos["second_read_out"]
                    ])

                exec_group.add_command(fix_cutadapt)

                # 4. Read data from first_read fifo
                with exec_group.add_pipeline() as fr_pigz_pipe:
                    # 4.1  command: Read from first_read fifos
                    cat = [self.get_tool('cat'),
                           temp_fifos["first_read_out"]]
                    # 4.2 Gzip output file
                    pigz = [self.get_tool('pigz'),
                            '--processes', str(self.get_cores()),
                            '--blocksize', self.get_option('pigz-blocksize'),
                            '--stdout']
                    # 4.3 command: Write to output file in 'dd-blocksize'
                    # chunks
                    fr_stdout_path = run.add_output_file(
                        "first_read",
                        "%s%s.fastq.gz" %
                        (run_id, read_types["first_read"]),
                        cc[run_id]["in/first_read"])
                    dd = [self.get_tool('dd'),
                          'obs=%s' % self.get_option('dd-blocksize'),
                          'of=%s' % fr_stdout_path]

                    fr_pigz_pipe.add_command(cat)
                    fr_pigz_pipe.add_command(pigz)
                    fr_pigz_pipe.add_command(dd)

                # 5. Read data from second_read fifo if there is one
                if cc.connection_exists("in/second_read"):
                    with exec_group.add_pipeline() as sr_pigz_pipe:
                        # 5.1  command: Read from first_read fifos
                        cat = [self.get_tool('cat'),
                               temp_fifos["second_read_out"]]
                        # 4.2 Gzip output file
                        pigz = [self.get_tool('pigz'),
                                '--processes',
                                str(self.get_cores()),
                                '--blocksize',
                                self.get_option('pigz-blocksize'),
                                '--stdout']
                        # 4.3 command: Write to output file in 'dd-blocksize'
                        # chunks
                        sr_stdout_path = run.add_output_file(
                            "second_read",
                            "%s%s.fastq.gz" %
                            (run_id, read_types["second_read"]),
                            cc[run_id]["in/second_read"])
                        dd = [self.get_tool('dd'),
                              'obs=%s' % self.get_option('dd-blocksize'),
                              'of=%s' % sr_stdout_path]

                        sr_pigz_pipe.add_command(cat)
                        sr_pigz_pipe.add_command(pigz)
                        sr_pigz_pipe.add_command(dd)
