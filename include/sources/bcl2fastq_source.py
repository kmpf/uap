import sys
from abstract_step import *
import os
import process_pool

class Bcl2FastqSource(AbstractStep):

    def __init__(self, pipeline):
        super(Bcl2FastqSource, self).__init__(pipeline)

        self.add_connection('out/configureBcl2Fastq_log_stderr')
        self.add_connection('out/make_log_stderr')
        self.add_connection('out/flow_cell_path')

        self.require_tool('configureBclToFastq.pl')
        self.require_tool('make')
        
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
        output_dir = '%s/Unaligned' % input_dir
        if 'output-dir' in options.keys():
            output_dir = options[option]
        # Define a new run
        run = self.new_run('read_demultiplexing')
        # Create new execution group for configureBclToFastq.pl
        bcl2Fastq_exec_group = run.new_exec_group()
        # Assemble complete command
        configureBcl2Fastq = ['configureBclToFastq.pl', '--input-dir',
                              '%s/Data/Intensities/BaseCalls/' % 
                              input_dir]

        for option in options.keys():
            if option in ['adapter-sequence', 'adapter-stringency',
                          'use-bases-mask', 'no-eamss',
                          'with-failed-reads', 'intensities-dir',
                          'positions-dir', 'positions-format',
                          'filter-dir', 'output-dir',
                          'sample-sheet', 'mismatches',
                          'fastq-cluster-count', 'ignore-missing-stats',
                          'ignore-missing-bcl', 'ignore-missing-control',
                          'tiles', 'flowcell-id'] :
                configureBcl2Fastq.extend([ '--%s' % option, options[option] ])
        # Create new command in execution group
        configureBcl2Fastq_command = bcl2Fastq_exec_group.new_command(
            configureBcl2Fastq,
            stderr_path = run.add_output_file(
                "configureBcl2Fastq_log_stderr",
                "bcl2fastq-log_stderr.txt", []) )
        print(configureBcl2Fastq)
        # Create new execution group for make
        make_exec_group = run.new_exec_group()
        make = ['make', '-C', '%s' % output_dir]
        make_command = make_exec_group.new_command(
            make,
            stderr_path = run.add_output_file(
                "make_log_stderr",
                "make-log_stderr.txt", []) )
        print(make)
        run.add_output_file("flow_cell_path",
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
