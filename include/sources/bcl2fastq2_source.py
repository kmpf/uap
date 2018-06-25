import sys
import logging
import os

from abstract_step import *

logger = logging.getLogger("uap_logger")

class Bcl2Fastq2Source(AbstractSourceStep):

    def __init__(self, pipeline):
        super(Bcl2Fastq2Source, self).__init__(pipeline)

        self.set_cores(10) # set # of cores for cluster, it is ignored if run locally

        #self.add_connection('out/configureBcl2Fastq_log_stderr')
        #self.add_connection('out/make_log_stderr')
        #self.add_connection('out/sample_sheet')
        self.add_connection('out/bcl2fastq2_log_stderr')

        self.require_tool('bcl2fastq')
        self.require_tool('mkdir')

        self.add_option('min-log-level', str, optional=True,
                        description='minimum log level recognized values: NONE, FATAL, '
                        'ERROR, WARNING, INFO, DEBUG, TRACE. (INFO)')
        self.add_option('runfolder-dir', str, optional=False,
                        description='path to runfolder directory (=./)')
        self.add_option('input-dir', str, optional=True,
                        description='path to input directory '
                        '(=<runfolder-dir>/Data/Intensities/BaseCalls/)')
        self.add_option('intensities-dir', str, optional=True,
                        description='path to intensities directory (=<input-dir>/../). '
                        'If intensities directory is specified, also input directory '
                        'must be specified.')
        self.add_option('output-dir', str, optional=False)
        self.add_option('interop-dir', str, optional=True,
                        description='path to demultiplexing statistics directory '
                        '(=<runfolder-dir>/InterOp/)')
        self.add_option('stats-dir', str, optional=True,
                        description='path to human-readable demultiplexing statistics '
                        'directory (=<runfolder-dir>/InterOp/)')
        self.add_option('reports-dir', str, optional=True,
                        description='path to reporting directory (=<output-dir>/Reports/)')
        self.add_option('sample-sheet', str, optional=True,
                        description='path to the sample sheet'
                        '(=<runfolder-dir>/SampleSheet.csv)')
        self.add_option('aggregated-tiles', str, optional=True,
                        description='tiles aggregation flag determining structure of '
                        'input files (=AUTO). recognized values: AUTO - Try to detect '
                        'correct setting. YES - Tiles are aggregated into single input '
                        'file. NO - There are separate input files for individual tiles')
        self.add_option('loading-threads', int, optional=True,
                        description='number of threads used for loading BCL data (=4)')
        self.add_option('processing-threads', int, optional=True,
                        description='number of threads used for processing demultipled '
                        'data (=100% of available CPUs)')
        self.add_option('writing-threads', int, optional=True,
                        description='number of threads used for writing FASTQ data '
                        '(=4)')
        self.add_option('tiles', str, optional=True,
                        description='Comma-separated list of regular expressions to select '
                        'only a subset of the tiles available in the flow-cell.Multiple '
                        'entries allowed, each applies to the corresponding base-calls.'
                        'For example: * to select all the tiles ending with 5 in all '
                        'lanes: tiles [0-9][0-9][0-9]5. * to select tile 2 in lane 1 '
                        'and all the tiles in the other lanes: tiles s_1_0002,s_[2-8]')
        self.add_option('minimum-trimmed-read-length', int, optional=True,
                        description='minimum read length after adapter trimming (=35)')
        self.add_option('use-bases-mask', str, optional=True,
                        description='Specifies how to use each cycle.')
        self.add_option('mask-short-adapter-reads', int, optional=True,
                        description='smallest number of remaining bases (after masking '
                        'bases below the minimum trimmed read length) below which whole '
                        'read is masked (=22)')
        self.add_option('adapter-stringency', float, optional=True,
                        description='adapter stringency (=0.9)')
        self.add_option('ignore-missing-bcls', bool, optional=True,
                        description='assume N for missing calls')
        self.add_option('ignore-missing-filter', bool, optional=True,
                        description='assume true for missing filters')
        self.add_option('ignore-missing-positions', bool, optional=True,
                        description='assume [0,i] for missing positions, where i is '
                        'incremented starting from 0')
        self.add_option('ignore-missing-controls', bool, optional=True,
                        description='assume 0 for missing controls')
        self.add_option('write-fastq-reverse-complement', bool, optional=True,
                        description='Generate FASTQs containing reverse complements of '
                        'actual data')
        self.add_option('with-failed-reads', bool, optional=True,
                        description='include non-PF clusters')
        self.add_option('create-fastq-for-index-reads', bool, optional=True,
                        description='create FASTQ files also for index reads')
        self.add_option('find-adapters-with-sliding-window', bool, optional=True,
                        description='find adapters with simple sliding window algorithm')
        self.add_option('no-bgzf-compression', bool, optional=True,
                        description='Turn off BGZF compression for FASTQ files')
        self.add_option('fastq-compression-level', int, optional=True,
                        description='Zlib compression level (1-9) used for FASTQ files (=4)')
        self.add_option('barcode-mismatches', bool, optional=True,
                        description='number of allowed mismatches per index multiple '
                        'entries, comma delimited entries, allowed; each entry is '
                        'applied to the corresponding index;last entry applies to all '
                        'remaining indices')
        self.add_option('no-lane-splitting', bool, optional=True,
                        description='Do not split fastq files by lane.')

    def runs(self, run_ids_connections_files):

        # Compile the list of options
        option_list = list()

        ## check all input folders for their existence
        path_options = ['input-dir', 'runfolder-dir', 'intensities-dir', 
                        'interop-dir', 'output-dir']
        set_path_options = [option for option in path_options if \
                            self.is_option_set_in_config(option)]

        for option in set_path_options:
            path = os.path.abspath(self.get_option(option))
            if os.path.isdir(path):
                option_list.append('--%s' % option)
                option_list.append(path)
            else:
                raise StandardError("No such directory for option '%s': %s" %
                                    (option, path))
        ## check if provided sample sheet exists
        file_option = 'sample-sheet'
        if self.is_option_set_in_config(file_option):
            file = os.path.abspath(self.get_option(file_option))
            if os.path.isfile(file):
                option_list.append('--%s' % file_option)
                option_list.append(file)
            else:
                raise StandardError("No such file for option '%s': %s" %
                                    (file_option, file))
        ### if the sample sheet is not provided in config, bcl2fastq automatically searches for
        ### it in the runfolder-dir, we need this file to split the process into runs per lane,
        ### which would speed up the whole demultiplexing process

        ## Get remaining options that are set in the configuration
        options = ["min-log-level", "stats-dir", "reports-dir", "aggregated-tiles",
                   "loading-threads", "demultiplexing-threads", "processing-threads",
                   "writing-threads", "tiles", "minimum-trimmed-read-length",
                   "use-bases-mask", "mask-short-adapter-reads", "adapter-stringency",
                   "ignore-missing-bcls", "ignore-missing-filter", "ignore-missing-positions",
                   "ignore-missing-controls", "write-fastq-reverse-complement",
                   "with-failed-reads", "create-fastq-for-index-reads",
                   "find-adapters-with-sliding-window", "no-bgzf-compression",
                   "fastq-compression-level", "barcode-mismatches", "no-lane-splitting"]
        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        for option in set_options:
            if self.get_option(option):
                option_list.append('--%s' % option)
            if not isinstance(self.get_option(option), bool):
                option_list.append(str(self.get_option(option)))

        # Declare a new run
        with self.declare_run('read_demultiplexing') as run:

            # the temporary output directory
            output_dir = run.get_output_directory_du_jour_placeholder()
            # Get path to Unaligned folder
            output_unaligned_dir = os.path.join(output_dir, 'Unaligned')

            if not os.path.exists(output_unaligned_dir):
                # Create new execution group for configureBclToFastq.pl
                with run.new_exec_group() as bcl2Fastq_exec_group:

                    bcl2fastq = [self.get_tool('bcl2fastq')]
                    bcl2fastq.extend(option_list)
                    bcl2fastq.extend(['--output-dir', output_unaligned_dir])

                    bcl2fastq_cmd = bcl2Fastq_exec_group.add_command(
                        bcl2fastq,
                        stderr_path = run.add_output_file("bcl2fastq2_log_stderr",
                                                          "bcl2fastq-log_stderr.txt", [])
                    )

                    logger.debug(" ".join(bcl2fastq))

                    # Create placeholder for Unaligned folder
#                    temp_output_dir = os.path.join(
#                        run.get_output_directory_du_jour_placeholder(),
#                        'Unaligned')
#
#                    # Assemble configureBclToFastq.pl command
#                    configureBcl2Fastq = [
#                        self.get_tool('configureBclToFastq.pl'),
#                        '--input-dir',
#                        os.path.join(input_dir, 'Data', 'Intensities',
#                                     'BaseCalls')
#                    ]
#
#                    configureBcl2Fastq.extend(option_list)
#
#                    configureBcl2Fastq.extend(
#                        ['--output-dir', temp_output_dir])
#
                    # Add command to execution group
#                    configureBcl2Fastq_command = bcl2Fastq_exec_group.add_command(
#                        configureBcl2Fastq,
#                        stderr_path = run.add_output_file(
#                            "configureBcl2Fastq_log_stderr",
#                            "bcl2fastq-log_stderr.txt", []) )
#                    logger.debug(" ".join(configureBcl2Fastq))
#                    # Create new execution group for make
#                    make_exec_group = run.new_exec_group()
#                    # Assemble make command
#                    make = [self.get_tool('make'), '-C', temp_output_dir, '-j', str(self.get_cores())]
#                    # Add make command to execution group
#                    make_exec_group.add_command(
#                        make,
#                        stderr_path = run.add_output_file(
#                            "make_log_stderr",
#                            "make-log_stderr.txt", []) )
#                    logger.debug(" ".join(make))
#                    logger.debug("Temporary output directory: %s" %
#                                 temp_output_dir)
#                    # Create new execution group to move Unaligned folder
#                    with run.new_exec_group() as mv_exec_group:
#                        # Assemble mkdir command, if output_dir does not exist
#
#                        mkdir = [self.get_tool('mkdir'), '-p',
#                                 output_unaligned_dir]
#                        logger.debug(" ".join(mkdir))
#                        # Add mkdir command to execution group
#
#                        mv_exec_group.add_command(mkdir)
#                        # Assemble mv command
#                        mv = [self.get_tool('mv'), temp_output_dir, output_dir]
#                        logger.debug(" ".join(mv))
#                        # Add mv command to execution group
#
#                        mv_command = mv_exec_group.add_command(mv)
#                        run.add_public_info("bcl2fastq-output-folder",
#                                            output_unaligned_dir)

