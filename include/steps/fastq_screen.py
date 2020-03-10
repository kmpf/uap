import sys
import logging
from abstract_step import AbstractStep
import os

logger = logging.getLogger('uap_logger')


class FastqScreen(AbstractStep):
    '''
    Fastq Screen (ver. 0.11.*)

    FastQ Screen allows you to screen a library of sequences in FastQ format
    against a set of sequence databases so you can see if the composition of
    the library matches with what you expect.

    '''

    def __init__(self, pipeline):
        super(FastqScreen, self).__init__(pipeline)
        cores = 10
        self.set_cores(cores)

        self.add_connection('in/first_read')
        self.add_connection('out/fqc_report')
        self.add_connection('out/fqc_image')
        self.add_connection('out/fqc_html')
        self.add_connection('out/tagged', optional=True)
        self.add_connection('out/tagged_filter', optional=True)
        self.add_connection('out/fastq_screen.conf', optional=True)
        self.add_connection('out/log_stdout')
        self.add_connection('out/log_stderr')

        self.add_option('config', str, optional=True,
                        description="Manually specify a location for the \
                        configuration.")

        self.add_option('databases', dict, optional=True,
                        description="Manually specify a location for the \
                        configuration. E.g.: fastq_screen: databases: Human: \
                        /path/to/human/bowtie2/index")

        self.add_option('keep config', bool, optional=True, default=False,
                        description="Keep the generated fastq_screen.conf \
                        for each run (only if databases are specified).")

        self.add_option('cores', int, default=cores)

        self.add_option('nohits', bool, default=False, optional=True,
                        description="Writes to a file the sequences that did \
                        not map to any of the specified genomes. This option \
                        is equivalent to specifying --tag --filter 0000 \
                        (number of zeros corresponds to the number of genomes \
                        screened).  By default the whole input file will be \
                        mapped, unless overridden by --subset.")

        self.add_option('subset', int, optional=True,
                        description="Don't use the whole sequence file, but \
                        create a temporary dataset of this specified number \
                        of reads. The dataset created will be of approximately \
                        (within a factor of 2) of this size. If the real \
                        dataset is smaller than twice the specified size \
                        then the whole dataset will be used. Subsets will \
                        be taken evenly from throughout the whole original \
                        dataset. By Default FastQ Screen runs with this \
                        parameter set to 100,000. To process an entire dataset \
                        however, adjust --subset to 0.")

        self.require_tool('fastq_screen')
        self.require_tool('bowtie2')
        self.require_tool('mv')
        self.require_tool('rm')
        self.require_tool('printf')

    def runs(self, cc):
        self.set_cores(self.get_option('cores'))

        if self.get_option('config') and self.get_option('databases'):
            raise StepError(self, "A config file and databases are specified.")

        if not self.get_option('config') and not self.get_option('databases'):
            raise StepError(self, "No config file or databases are specified.")

        if self.get_option('config'):
            logger.warning(
                '[%s] Using a config file is deprecated. '
                'Please specify databases instead.' %
                self.get_step_name())
            config_file = os.path.abspath(self.get_option('config'))
        else:
            conf_data = [
                'BOWTIE2 %s' % self.get_tool('bowtie2'),
                'THREADS %d' % self.get_cores()
            ]
            for db in sorted(self.get_option('databases').items()):
                conf_data.append('DATABASE %s %s BOWTIE2' % db)

        for run_id in cc.keys():
            run = self.declare_run(run_id)
            if not self.get_option('config'):
                if self.get_option('keep config'):
                    config_file = run.add_output_file('fastq_screen.conf',
                                                      'fastq_screen.conf', [])
                else:
                    config_file = run.add_temporary_file('fastq_screen.conf')
                write_conf = [self.get_tool('printf'), '\n'.join(conf_data)]
                execg = run.new_exec_group()
                execg.add_command(write_conf, stdout_path=config_file)
            for input_path in cc[run_id]['in/first_read']:
                file_name = os.path.basename(input_path).rstrip(".fastq.gz")
                # prepare output files
                file_pattern = "%s_screen.txt" % (file_name)
                run.add_output_file("fqc_report", file_pattern, [input_path])

                file_pattern = "%s_screen.png" % (file_name)
                run.add_output_file("fqc_image", file_pattern, [input_path])

                file_pattern = "%s_screen.html" % (file_name)
                run.add_output_file("fqc_html", file_pattern, [input_path])

                file_pattern = "%s-fastqscreen-log_stdout.txt" % (file_name)
                log_stdout = run.add_output_file("log_stdout",
                                                 file_pattern,
                                                 [input_path])

                file_pattern = "%s-fastqscreen-log_stderr.txt" % (file_name)
                log_stderr = run.add_output_file("log_stderr",
                                                 file_pattern,
                                                 [input_path])

                # build fastq_screen command
                fastq_screen_exec_group = run.new_exec_group()
                fastq_screen = [self.get_tool('fastq_screen'),
                                '-conf', config_file]

                if self.get_option('subset'):
                    fastq_screen.extend(['--subset',
                                         str(self.get_option('subset'))])

                if self.get_option('nohits'):
                    file_pattern = "%s.tagged.fastq.gz" % (file_name)
                    run.add_output_file("tagged", file_pattern, [input_path])

                    file_pattern = "%s.tagged_filter.fastq.gz" % (file_name)
                    run.add_output_file("tagged_filter", file_pattern,
                                        [input_path])

                    fastq_screen.extend(['--nohits'])

                fastq_screen.extend(['--outdir', '.', input_path])
                fastq_screen_exec_group.add_command(
                    fastq_screen, stdout_path=log_stdout, stderr_path=log_stderr)
