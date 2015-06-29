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
        self.add_connection('out/unaligned_folder_path')

        self.require_tool('configureBclToFastq.pl')
        self.require_tool('make')
        self.require_tool('mv')
        self.require_tool('mkdir')
        
#        self.add_option('paired_end', bool)
#        self.add_option('first_read', str, default = "_R1",
#            description = "Part of the file name that marks all "
#                "files containing sequencing data of the first read. "
#                "Example: '_R1.fastq' or '_1.fastq'")
#
#        self.add_option('second_read', str, default = "_R2",
#            description = "Part of the file name that marks all "
#                "files containing sequencing data of the second read. "
#                "Example: 'R2.fastq' or '_2.fastq'")
#

        self.add_option('input-dir', str, optional=False, description="file URL")
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
        self.add_option('output-dir', str, optional=True)
        self.add_option('sample-sheet', str, optional=True)
        self.add_option('mismatches', str, optional=True)
        self.add_option('fastq-cluster-count', str, optional=True)
        self.add_option('ignore-missing-stats', str, optional=True)
        self.add_option('ignore-missing-bcl', str, optional=True)
        self.add_option('ignore-missing-control', str, optional=True)
        self.add_option('tiles', str, optional=True)
        self.add_option('flowcell-id', str, optional=True)


    def runs(self, run_ids_connections_files):
        '''
        
        '''
        options = self.get_options()
        input_dir = self.get_option('input-dir')
        # Create path to Unaligned folder ...
        output_dir = '%s%s' % (input_dir, 'Unaligned')
        # ... or use the given path
        if 'output-dir' in options.keys():
            output_dir = options['output-dir']

        # Create placeholder for Unaligned folder
        temp_output_dir = os.path.join(
            self.get_output_directory_du_jour_placeholder(),
            'Unaligned')
        # Define a new run
        run = self.new_run('read_demultiplexing')
        # Create new execution group for configureBclToFastq.pl
        bcl2Fastq_exec_group = run.new_exec_group()
        # Assemble configureBclToFastq.pl command
        configureBcl2Fastq = [self.get_tool('configureBclToFastq.pl'),
                              '--input-dir',
                              '%s/Data/Intensities/BaseCalls/' % 
                              input_dir]

        for option in options.keys():
            if option in ['adapter-sequence', 'adapter-stringency',
                          'use-bases-mask', 'no-eamss',
                          'with-failed-reads', 'intensities-dir',
                          'positions-dir', 'positions-format',
                          'filter-dir', 'sample-sheet', 'mismatches',
                          'fastq-cluster-count', 'ignore-missing-stats',
                          'ignore-missing-bcl', 'ignore-missing-control',
                          'tiles', 'flowcell-id'] :
                configureBcl2Fastq.extend([ '--%s' % option, options[option] ])
        configureBcl2Fastq.extend(
                    ['--output-dir', temp_output_dir])
        # Add command to execution group
        configureBcl2Fastq_command = bcl2Fastq_exec_group.new_command(
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
        if not os.path.exists(output_dir):
            mkdir = [self.get_tool('mkdir'), output_dir]
            logger.debug(" ".join(mkdir))
            # Add mkdir command to execution group
            mv_exec_group.new_command(mkdir)
        # Assemble mv command
        mv = [self.get_tool('mv')]
        mv.extend([temp_output_dir, output_dir])
        logger.debug(" ".join(mv))
        # Add mv command to execution group
        mv_command = mv_exec_group.new_command(mv)
        # Add output file
        run.add_output_file("unaligned_folder_path",
                            output_dir, [])

#        baseDir=/data/bioinf/Data/141121_SN928_0101_AC39TWACXX/
#        outDir=/data/bioinf/Data/141121_SN928_0101_AC39TWACXX_Keep/\
#            Unaligned_CentOS6/
#        sampleSheet=./SampleSheet2_corrected.csv
#        configureBclToFastq.pl --use-bases-mask Y101,I6n,Y51n --mismatches 1 \
#            --input-dir ${baseDir}/Data/Intensities/BaseCalls/ --sample-sheet \
#            ${sampleSheet} --fastq-cluster-count 1000000000 --output-dir \
#            ${outDir}
#        make -j 20 -C Unaligned_CentOS6
