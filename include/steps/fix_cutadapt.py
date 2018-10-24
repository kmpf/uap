import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class FixCutadapt(AbstractStep):
    '''
    This step takes FASTQ data and removes both reads of a paired-end read, if
    one of them has been completely removed by cutadapt (or any other software).
    '''
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "2M")
        self.add_option('pigz-blocksize', str, optional = True, default = "2048")

        # Step was tested for cat (GNU coreutils) release 8.25
        self.require_tool('cat')
        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        self.require_tool('fix_cutadapt')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

    def runs(self, run_ids_connections_files):

        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                temp_fifos = dict()
                exec_group = run.new_exec_group()
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]
                    temp_fifos["%s_in" % read] = None
                    temp_fifos["%s_out" % read] = None
                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)

                    elif len(input_paths) != 1:
                        logger.error("Expected single input file. Found files "
                                     "%s for run: %s" % (input_paths, run_id))
                        sys.exit(1)
                    else:
                        # 1. Create temporary fifos
                        # 1.1 Input fifo
                        temp_fifos["%s_in" % read] = run.add_temporary_file(
                            "in-fifo-%s" %
                            os.path.basename(input_paths[0]) )
                        mkfifo_in = [self.get_tool('mkfifo'),
                                     temp_fifos["%s_in" % read]]
                        exec_group.add_command(mkfifo_in)
                        # 1.2 Output fifo
                        temp_fifos["%s_out" % read] = run.add_temporary_file(
                            "%s-out-fifo" % read)
                        mkfifo_out = [self.get_tool('mkfifo'),
                                      temp_fifos["%s_out" % read]]
                        exec_group.add_command(mkfifo_out)

                        # 2. Output files to fifo
                        if input_paths[0].endswith('fastq.gz'):
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in 'dd-blocksize' chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'ibs=%s' % self.get_option('dd-blocksize'),
                                    'if=%s' % input_paths[0]
                                ]
                                # 2.2 command: Uncompress file to fifo
                                pigz = [self.get_tool('pigz'),
                                        '--processes', str(self.get_cores()),
                                        '--decompress',
                                        '--blocksize', self.get_option('pigz-blocksize'),
                                        '--stdout']
                                # 2.3 Write file in 'dd-blocksize' chunks to fifo
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
                            #              write to fifo in 'dd-blocksize' chunks
                            dd_in = [
                                self.get_tool('dd'),
                                'bs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_paths[0],
                                'of=%s' % temp_fifos["%s_in" % read]
                            ]
                            exec_group.add_command(dd_in)
                        else:
                            logger.error("File %s does not end with any "
                                         "expected suffix (fastq.gz or fastq). "
                                         "Please fix that issue." % input_path)
                            sys.exit(1)
                # 3. Start fix_cutadapt
                fix_cutadapt = [self.get_tool('fix_cutadapt'),
                                temp_fifos["first_read_in"], 
                                temp_fifos["first_read_out"] ]
                if temp_fifos["second_read_in"] != None and \
                   temp_fifos["second_read_out"] != None:
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
                    # 4.3 command: Write to output file in 'dd-blocksize' chunks
                    fr_stdout_path = run.add_output_file(
                        "first_read",
                        "%s%s.fastq.gz" %
                        (run_id, read_types["first_read"]),
                        run_ids_connections_files[run_id]["in/first_read"])
                    dd = [self.get_tool('dd'),
                          'obs=%s' % self.get_option('dd-blocksize'),
                          'of=%s' % fr_stdout_path]

                    fr_pigz_pipe.add_command(cat)
                    fr_pigz_pipe.add_command(pigz)
                    fr_pigz_pipe.add_command(dd)

                # 5. Read data from second_read fifo if there is one
                if temp_fifos["second_read_in"] != None and \
                   temp_fifos["second_read_out"] != None:
                    with exec_group.add_pipeline() as sr_pigz_pipe:
                        # 5.1  command: Read from first_read fifos
                        cat = [self.get_tool('cat'),
                               temp_fifos["second_read_out"]]
                        # 4.2 Gzip output file
                        pigz = [self.get_tool('pigz'),
                                '--processes', str(self.get_cores()),
                                '--blocksize', self.get_option('pigz-blocksize'),
                                '--stdout']
                        # 4.3 command: Write to output file in 'dd-blocksize' chunks
                        sr_stdout_path = run.add_output_file(
                            "second_read",
                            "%s%s.fastq.gz" %
                            (run_id, read_types["second_read"]),
                            run_ids_connections_files[run_id]["in/second_read"])
                        dd = [self.get_tool('dd'),
                              'obs=%s' % self.get_option('dd-blocksize'),
                              'of=%s' % sr_stdout_path]

                        sr_pigz_pipe.add_command(cat)
                        sr_pigz_pipe.add_command(pigz)
                        sr_pigz_pipe.add_command(dd)
