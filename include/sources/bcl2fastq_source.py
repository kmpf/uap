import sys
import logging
import os

from abstract_step import *

logger = logging.getLogger("uap_logger")

class Bcl2FastqSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(Bcl2FastqSource, self).__init__(pipeline)

        self.add_connection('out/configureBcl2Fastq_log_stderr')
        self.add_connection('out/make_log_stderr')
        self.add_connection('out/sample_sheet')

        self.require_tool('configureBclToFastq.pl')
        self.require_tool('make')
        self.require_tool('mv')
        self.require_tool('mkdir')

        self.add_option('input-dir', str, optional=False, description="file URL")
        self.add_option('output-dir', str, optional=True)
        self.add_option('adapter-sequence', str, optional=True)
        self.add_option('adapter-stringency', str, optional=True)
        self.add_option('use-bases-mask', str, optional=True, description=
                        'Conversion mask characters:'
                        '- Y or y: use'
                        '- N or n: discard'
                        '- I or i: use for indexing'
                        'If not given, the mask will be guessed from the'
                        'RunInfo.xml file in the run folder.'
                        'For instance, in a 2x76 indexed paired end run, the'
                        'mask *Y76,I6n,y75n* means: "use all 76 bases from the'
                        'first end, discard the last base of the indexing read,'
                        'and use only the first 75 bases of the second end".')
        self.add_option('no-eamss', str, optional=True)
        self.add_option('with-failed-reads', str, optional=True)
        self.add_option('intensities-dir', str, optional=True)
        self.add_option('positions-dir', str, optional=True)
        self.add_option('positions-format', str, optional=True)
        self.add_option('filter-dir', str, optional=True)
        self.add_option('sample-sheet', str, optional=True)
        self.add_option('mismatches', int, optional=True)
        self.add_option('fastq-cluster-count', int, optional=True)
        self.add_option('ignore-missing-stats', str, optional=True)
        self.add_option('ignore-missing-bcl', str, optional=True)
        self.add_option('ignore-missing-control', str, optional=True)
        self.add_option('tiles', str, optional=True)
        self.add_option('flowcell-id', str, optional=True)


    def runs(self, run_ids_connections_files):
        '''
        
        '''
        input_dir = self.get_option('input-dir')
        # Create path to Unaligned folder ...
        output_dir = os.path.join(input_dir)
        # ... or use the given path
        if self.is_option_set_in_config('output-dir'):
            output_dir = self.get_option('output-dir')

        # Create placeholder for Unaligned folder
        temp_output_dir = os.path.join(
            self.get_output_directory_du_jour_placeholder(),
            'Unaligned')

        # Declare a new run
        with self.declare_run('read_demultiplexing') as run:
            # Create new execution group for configureBclToFastq.pl
            bcl2Fastq_exec_group = run.new_exec_group()
            # Assemble configureBclToFastq.pl command
            configureBcl2Fastq = [self.get_tool('configureBclToFastq.pl'),
                                  '--input-dir',
                                  os.path.join(input_dir, 'Data', 'Intensities', 
                                               'BaseCalls')]

            options = ['adapter-sequence', 'adapter-stringency', 'use-bases-mask',
                       'no-eamss', 'with-failed-reads', 'intensities-dir',
                       'positions-dir', 'positions-format', 'filter-dir',
                       'sample-sheet', 'mismatches', 'fastq-cluster-count',
                       'ignore-missing-stats', 'ignore-missing-bcl', 
                       'ignore-missing-control', 'tiles', 'flowcell-id']
            set_options = [option for option in options if \
                           self.is_option_set_in_config(option)]
            for option in set_options:
                configureBcl2Fastq.extend([ '--%s' % option, 
                                            str(self.get_option(option)) ])
            configureBcl2Fastq.extend(
                ['--output-dir', temp_output_dir])
            # Add command to execution group
            configureBcl2Fastq_command = bcl2Fastq_exec_group.add_command(
                configureBcl2Fastq,
                stderr_path = run.add_output_file(
                    "configureBcl2Fastq_log_stderr",
                    "bcl2fastq-log_stderr.txt", []) )
            logger.debug(" ".join(configureBcl2Fastq))
            # Create new execution group for make
            make_exec_group = run.new_exec_group()
            # Assemble make command
            make = [self.get_tool('make'), '-C', temp_output_dir]
            # Add make command to execution group
            make_exec_group.new_command(
                make,
                stderr_path = run.add_output_file(
                    "make_log_stderr",
                    "make-log_stderr.txt", []) )
            logger.debug(" ".join(make))
            logger.debug("Temporary output directory: %s" % temp_output_dir)
            # Create new execution group to move Unaligned folder
            mv_exec_group = run.new_exec_group()
            # Assemble mkdir command, if output_dir does not exist
            output_unaligned_dir = os.path.join(output_dir, 'Unaligned')
            if not os.path.exists(output_unaligned_dir):
                mkdir = [self.get_tool('mkdir'), '-p', output_unaligned_dir]
                logger.debug(" ".join(mkdir))
                # Add mkdir command to execution group

                mv_exec_group.add_command(mkdir)
                # Assemble mv command
                mv = [self.get_tool('mv', temp_output_dir, output_dir])
                logger.debug(" ".join(mv))
                # Add mv command to execution group

                mv_command = mv_exec_group.add_command(mv)
                run.add_public_info("bcl2fastq-output-folder",
                                    output_unaligned_dir)
            else:
                logger.warning("Directory: %s already exists!" %
                               output_unaligned_dir)
                logger.warning("If you WANT to restart the analysis either "
                               "remove the directory or choose a new output "
                               "directory. ")
                logger.warning("If you DON'T want to restart the analysis "
                               "everything is fine.")

