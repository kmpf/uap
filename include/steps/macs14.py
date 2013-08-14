import sys
from abstract_step import *
import pipeline
import re
import process_pool
import yaml

class Macs14(AbstractStep):
    
    def __init__(self, pipeline):
        super(Macs14, self).__init__(pipeline)
        
        self.set_cores(4)

        self.add_connection('in/alignments')
#        self.add_connection('out/peaks')
        self.add_connection('out/log')
 
        self.require_tool('macs14')
        self.require_tool('cat4m')
        self.require_tool('pigz')

        self.add_option('control', list , default=list())
        self.add_option('format', str, default='AUTO',
                        choices=['ELAND', 'ELANDMULTI', 'ELANDMULTIPET',
                                 'ELANDEXPORT', 'SAM', 'BAM', 'BOWTIE'])
        self.add_option('genome_size', str, default='hs')


    def declare_runs(self):
        
        control_files = list()
        validated_control_samples = list()
        control_samples = list()
        if self.get_option('control'):
            control_samples.extend(self.get_option('control'))

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
#                print(input_paths)
                run.add_output_file('log', '%s-macs14-log.txt' % run_id, input_paths)
                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                run.add_private_info('treatment_files', input_paths)

                if run_id in control_samples:
                    validated_control_samples.append(run_id)
                    control_files.extend(input_paths)

        for run_id in self.get_run_ids():
            run = self.get_run(run_id)
            if run != None:
                if validated_control_samples and control_files:
#                    print("Control files for runID %s: %s" % (run_id, control_files))
                    run.add_private_info('control_files', control_files)
#                    print("Control samples for runID %s: %s" % (run_id, validated_control_samples))
                    run.add_public_info('control_samples', validated_control_samples)
                
    
    def execute(self, run_id, run):

        with process_pool.ProcessPool(self) as pool:

            with pool.Pipeline(pool) as pipeline:

                macs14 = [self.get_tool('macs14'), '--treatment']
                if not run.has_private_info('treatment_files'):
                    raise StandardError("No treatment files for %s to analyse with macs14" % run_id)
                macs14.extend(run.get_private_info('treatment_files'))

                # if we do have control data use it
                if run.has_private_info('control_files'):
                    macs14.extend(['--control'])
                    macs14.extend(run.get_private_info('control_files'))

                macs14.extend([
                        '--format', self.get_option('format'),
                        '--name', run_id
                        ])
                pipeline.append(macs14)

                # MACS14 use without control sample
                # macs14 --treatment=$j --format=BED --name $(basename $j .bed) -S >
                #        $(basename $j .bed).log

                # MACS14 use with control sample
                #     macs14 --treatment=$i --control=../rabbit_IgG_TGACCA_L001_001.bed \
                #            --format=BED --name=$(basename $i .bed) -S > $(basename $i .bed).log
