import sys
from abstract_step import *
import misc
import process_pool
import yaml


class Segemehl(AbstractStep):


    def __init__(self, pipeline):
        super(Segemehl, self).__init__(pipeline)

        self.set_cores(12)
        
        self.add_connection('in/reads')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/log')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('segemehl')

        self.add_option('genome', str)
        self.add_option('index', str)
        

    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            
            with self.declare_run(run_id) as run:
                read_files = misc.assign_strings(input_paths, ['R1','R2'])
                run.add_output_file('alignments', '%s-segemehl-results.sam.gz' % run_id, input_paths)
                run.add_output_file('unmapped', '%s-segemehl-unmapped.fastq.gz' % run_id, input_paths)
                run.add_output_file('log', '%s-segemehl-log.txt' % run_id, input_paths)
                
                run.add_private_info('R1-in', read_files['R1'])
                run.add_private_info('R2-in', read_files['R2'])


    def execute(self, run_id, run):

        with process_pool.ProcessPool(self) as pool:
            
            fifo_path_genome = pool.get_temporary_fifo('segemehl-genome-fifo', 'input')
            fifo_path_unmapped = pool.get_temporary_fifo('segemehl-unmapped-fifo', 'output')
            
            pool.launch([self.get_tool('cat4m'), self.get_option('genome'), '-o', fifo_path_genome])
            
            with pool.Pipeline(pool) as pipeline:
            
                q = run.get_private_info('R1-in')
                p = run.get_private_info('R2-in')
                    
                segemehl = [
                    self.get_tool('segemehl'),
                    '-d', fifo_path_genome,
                    '-i', self.get_option('index'),
                    '-q', q,
                    '-p', p,
                    '-u', fifo_path_unmapped,
                    '-H', '1',
                    '-t', '11',
                    '-s', '-S',
                    '-D', '0',
                    '-o', '/dev/stdout'
                ]
                
                pigz = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
                
                pipeline.append(segemehl, stderr_path = run.get_single_output_file_for_annotation('log'))
                pipeline.append(pigz, stdout_path = run.get_single_output_file_for_annotation('alignments'))
                
            with pool.Pipeline(pool) as pipeline:
                pigz = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
                
                pipeline.append([self.get_tool('cat4m'), fifo_path_unmapped])
                pipeline.append(pigz, stdout_path = run.get_single_output_file_for_annotation('unmapped'))
                
