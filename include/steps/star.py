import sys
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Star(AbstractStep):
    '''

    '''
    def __init__(self, pipeline):
        super(Star, self).__init__(pipeline)

        # input connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        #self.add_connection('in/index')

        # output connections
        self.add_connection('out/aligned')
        self.add_connection('out/log.final')
        self.add_connection('out/log.out')
        self.add_connection('out/log.progess')
        self.add_connection('out/sj.out')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        # required tools
        self.require_tool('star')
        self.require_tool('rm')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                                    engine and threads ie")

        self.add_option('genomeDir', str, optional=True, default=None,
                        description="path to the directory where genome files \
                                are stored (if runMode!=generateGenome) or will be \
                                generated (if runMode==generateGenome)")

        self.add_option('readFilesCommand', str, optional=True, default=None,
                        description="command line to execute for each of the \
                        input file. This command should generate FASTA or \
                        FASTQ text and send it to stdout. \
                        For example: zcat - to uncompress .gz files, \
                        bzcat - to uncompress .bz2 files, etc.")

        self.add_option('runThreadN', int, optional=True, default=1,
                        description="number of threads to run STAR")

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

                star = [self.get_tool('star')]

                # get genomeDir from config or from input files
                if self.is_option_set_in_config('genomeDir'):
                    genome_dir = str(self.get_option('genomeDir'))
                else:
                    if 'in/genome_dir' not in run_ids_connections_files[run_id]:
                        logger.error('Required parameter "GenomDir" wasnt found!')
                        StandardError()
                    genome_dir = run_ids_connections_files[run_id]['in/genome_dir'][0]

                star.extend(['--genomeDir', genome_dir])

                out_path = run.get_output_directory_du_jour_placeholder()
                star.extend(['--outFileNamePrefix', out_path + '/'])

                if self.is_option_set_in_config('readFilesCommand'):
                    star.extend(['--readFilesCommand', self.get_option('readFilesCommand')])

                if self.is_option_set_in_config('cores'):
                    star.extend(['--runThreadN', str(self.get_option('runThreadN'))])

                star.append('--readFilesIn')
                star.extend(input_fileset)

                stderr_file = "%s-star-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, input_fileset)
                stdout_file = "%s-star-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, input_fileset)

                run.add_output_file("aligned", "Aligned.out.sam", input_fileset)
                run.add_output_file("log.final", "Log.final.out", input_fileset)
                run.add_output_file("log.out", "Log.out", input_fileset)
                run.add_output_file("log.progess", "Log.progress.out", input_fileset)
                run.add_output_file("sj.out", "SJ.out.tab", input_fileset)

                star_eg = run.new_exec_group()
                star_eg.add_command(star, stdout_path=log_stdout,
                                    stderr_path=log_stderr)

                # delete _STARtmp @ tmp directory if its not removed by STAR
                # ADDITIONAL COMMENT: this lines makes the result ambiguous 
                # because of the tmp_dir which uses the current time so
                # the uap want to produces the results repeatedly
                #tmp_dir = run.get_temp_output_directory() + '/_STARtmp'
                #eg = run.new_exec_group()
                #eg.add_command([self.get_tool('rm'), '-r', tmp_dir])
