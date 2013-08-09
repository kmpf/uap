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


        self.add_option('control', str, default='')
        self.add_option('format', str, default='AUTO')
        self.add_option('genome_size', str, default='hs')

    def declare_runs(self):
        
        control_file = str

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                print(input_paths)
                run.add_output_file('log', '%s-macs14-log.txt' % run_id, input_paths)

                if self.get_option('control'):
                    control = self.get_option('control')
                    run.add_public_info('control_sample', control)
                    print(self.get_option('control'))
                    
                    if control in run_id and not run.has_public_info('control_sample'):
                        run.add_public_info('control_sample', control)
                        if len(input_paths) != 1:
                            raise StandardError("Control group %s consists of zero or multiple files %s" % control, input_paths.join(', '))
                        control_file = input_paths[0]

        for run_id in self.get_run_ids():
            run = self.get_run(run_id)
            if run != None:
                if control_file:
                    run.add_private_info('control_file', control_file)
                
    
    def execute(self, run_id, run):

        with process_pool.ProcessPool(self) as pool:
            fifo_path_control = None
            fifo_path_sample = None
            if run.has_prublic_info('control_sample'):
                with pool.Pipeline(pool) as pipeline:
                    fifo_path_control = pool.get_temporary_fifo('macs14-control-fifo', 'input')
                    cat4m = [self.get_tool('cat4m'), run.get_private_info('control_file')]
                    pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1',
                            fifo_path_control]
                    pipeline.append(cat4m)
                    pipeline.append(pigz)

            with pool.Pipeline(pool) as pipeline:
                sample_file = run.get_single_input_file_for_annotation('alignments')
                fifo_path_sample = pool.get_temporary_fifo('%s-fifo' % sample_file, 'input')
                cat4m = [self.get_tool('cat4m'), sample_file]
                pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1',
                        fifo_path_sample]
                pipeline.append(cat4m)
                pipeline.append(pigz)

            with pool.Pipeline(pool) as pipeline:

                macs14 = [
                    self.get_tool('macs14'),
                    '--treatment', fifo_path_sample
                    ]

                # if we do have control data use it
                if fifo_path_control != None:
                    macs14.extend(['--control', fifo_path_control])
                
                macs14.extend([
                        '--format', self.get_option('format'),
                        '--name', run.get_private_info()
                        ])
                pipeline.append(macs14, stdout_path=run.get_private_info())

                # MACS14 use without control sampel
                # macs14 --treatment=$j --format=BED --name $(basename $j .bed) -S >
                #        $(basename $j .bed).log

                # MACS14 use with control sample
                #     macs14 --treatment=$i --control=../rabbit_IgG_TGACCA_L001_001.bed \
                #            --format=BED --name=$(basename $i .bed) -S > $(basename $i .bed).log
