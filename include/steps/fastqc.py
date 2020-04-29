import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Fastqc(AbstractStep):
    '''
    The fastqc step  is a wrapper for the fastqc tool. It generates some quality
    metrics for fastq files. For this specific instance only the zip archive is
    preserved.

    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

    Tested fastqc release: 0.11.2
    '''

    def __init__(self, pipeline):
        super(Fastqc, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read', optional=True)
        self.add_connection('out/first_read_fastqc_report')
        self.add_connection('out/first_read_fastqc_report_webpage')
        self.add_connection('out/first_read_log_stderr')
        self.add_connection('out/second_read_fastqc_report', optional=True)
        self.add_connection(
            'out/second_read_fastqc_report_webpage',
            optional=True)
        self.add_connection('out/second_read_log_stderr', optional=True)
        self.require_tool('fastqc')
        # Step was tested for mkdir (GNU coreutils) release 8.25
        self.require_tool('mkdir')
        # Step was tested for mv (GNU coreutils) release 8.25
        self.require_tool('mv')

        # [Options for 'fastqc':]
        # -o|--outdir is set by uap
        # --extract and --noextract are controlled by uap since it influences the output
        self.add_option(
            'casava',
            bool,
            optional=True,
            description="Files come from raw casava output. Files in the same "
            "sample group (differing only by the group number) will be analysed "
            "as a set rather than individually. Sequences with the filter flag "
            "set in the header will be excluded from the analysis. Files must "
            "have the same names given to them by casava (including being gzipped "
            "and ending with .gz) otherwise they won't be grouped together correctly.")
        self.add_option('nofilter', bool, optional=True,  # default=False,
                        description="If running with --casava then do not remove read flagged by "
                        "casava as poor quality when performing the QC analysis.")
        self.add_option(
            'java',
            str,
            optional=True,
            description="Provides the full path to the java binary you want to use "
            "to launch fastqc. If not supplied then java is assumed to be in your path.")
        self.add_option(
            'nogroup',
            bool,
            optional=True,
            description="Disable grouping of bases for reads >50bp. All reports will "
            "show data for every base in the read.  WARNING: Using this option will "
            "cause fastqc to crash and burn if you use it on really long reads, and "
            "your plots may end up a ridiculous size. You have been warned!")
        self.add_option(
            'format',
            str,
            optional=True,
            choices=[
                'bam',
                'sam',
                'bam_mapped',
                'sam_mapped',
                'fastq'],
            description="Bypasses the normal sequence file format detection and forces "
            "the program to use the specified format.  Valid formats are bam,sam, "
            "bam_mapped,sam_mapped and fastq")
        self.add_option(
            'threads',
            int,
            optional=True,
            description="Specifies the number of files which can be processed "
            "simultaneously. Each thread will be allocated 250MB of memory so "
            "you should not run more threads than your available memory will cope "
            "with, and not more than 6 threads on a 32 bit machine")
        self.add_option(
            'contaminants',
            str,
            optional=True,
            description="Specifies a non-default file which contains the list of "
            "contaminants to screen overrepresented sequences against. The file must "
            "contain sets of named contaminants in the form name[tab]sequence.  Lines "
            "prefixed with a hash will be ignored.")
        self.add_option(
            'adapters',
            str,
            optional=True,
            description="Specifies a non-default file which contains the list of "
            "adapter sequences which will be explicity searched against the library. "
            "The file must contain sets of named adapters in the form name[tab]sequence. "
            "Lines prefixed with a hash will be ignored.")
        self.add_option(
            'limits',
            str,
            optional=True,
            description="Specifies a non-default file which contains a set of criteria "
            "which will be used to determine the warn/error limits for the various "
            "modules.  This file can also be used to selectively remove some modules "
            "from the output all together.  The format needs to mirror the default "
            "limits.txt file found in the Configuration folder.")
        self.add_option(
            'kmers',
            int,
            optional=True,
            description="Specifies the length of Kmer to look for in the Kmer content "
            "module. Specified Kmer length must be between 2 and 10. Default length is "
            "7 if not specified.")
        self.add_option(
            'dir',
            str,
            optional=True,
            description="Selects a directory to be used for temporary files written "
            "when generating report images. Defaults to system temp directory if not "
            "specified.")

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional=True, default="2M")
        self.add_option('pigz-blocksize', str, optional=True, default="2048")

    def runs(self, run_ids_connections_files):

        options = [
            'casava',
            'nofilter',
            'java',
            'nogroup',
            'format',
            'contaminants',
            'adapters',
            'limits',
            'kmers',
            'dir']

        set_options = [option for option in options if
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        if 'threads' not in set_options:
            option_list.append('--threads')
            option_list.append(str(self.get_cores()))
        else:
            self.set_cores(self.get_option('threads'))

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id].get(
                        connection)
                    if input_paths:
                        for input_path in input_paths:
                            # Get base name of input file
                            root, ext = os.path.splitext(os.path.basename(
                                input_path))
                            if os.path.basename(input_path).endswith(
                                    ('.fq.gz', '.fq.gzip', '.fastq.gz',
                                     '.fastq.gzip')):
                                parts = os.path.basename(input_path).split('.')
                                root = '.'.join(parts[:-2])
                                ext = '.'.join(parts[-2:])

                            # Create temporary output directory
                            temp_dir = run.add_temporary_directory(
                                "%s" % root)
                            mkdir_exec_group = run.new_exec_group()
                            mkdir = [self.get_tool('mkdir'), temp_dir]
                            mkdir_exec_group.add_command(mkdir)
                            # 1. Run fastqc for input file
                            fastqc_exec_group = run.new_exec_group()
                            fastqc = [self.get_tool('fastqc'), '--noextract']
                            fastqc.extend(option_list)
                            fastqc.extend(['-o', temp_dir])
                            fastqc.append(input_path)
                            fastqc_command = fastqc_exec_group.add_command(
                                fastqc,
                                stderr_path=run.add_output_file(
                                    "%s_log_stderr" % read,
                                    "%s%s-fastqc-log_stderr.txt" %
                                    (run_id, read_types[read]),
                                    [input_path]))
                            # 2. Move fastqc results to final destination
                            mv_exec_group = run.new_exec_group()
                            mv1 = [self.get_tool('mv'),
                                   os.path.join(temp_dir,
                                                ''.join([root,
                                                         '_fastqc.zip'])),
                                   run.add_output_file(
                                "%s_fastqc_report" % read,
                                "%s%s-fastqc.zip" %
                                (run_id, read_types[read]),
                                [input_path])]
                            mv2 = [self.get_tool('mv'),
                                   os.path.join(temp_dir,
                                                ''.join([root,
                                                         '_fastqc.html'])),
                                   run.add_output_file(
                                "%s_fastqc_report_webpage" % read,
                                "%s%s-fastqc.html" %
                                (run_id, read_types[read]),
                                [input_path])]

                            mv_exec_group.add_command(mv1)
                            mv_exec_group.add_command(mv2)
