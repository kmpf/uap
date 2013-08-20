import sys
from abstract_step import *
import misc
import process_pool
import yaml

class RandomSampleFastq(AbstractStep):
    
    def __init__(self, pipeline):
        super(RandomSampleFastq, self).__init__(pipeline)
                        
        self.set_cores(3)
        #set mem  #necessary
        self.add_connection('in/reads')
        self.add_connection('out/reads')

        self.require_tool('pigz')
        self.require_tool('random_sample_fastq')

        self.add_option('sample_size', int, default = 1000)
        
    def declare_runs(self):

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            with self.declare_run(run_id) as run:
                is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
                if is_paired_end:
                    run.add_output_file('reads', "%s-subsample-R1.fastq.gz" % run_id, input_paths)
                    run.add_output_file('reads', "%s-subsample-R2.fastq.gz" % run_id, input_paths)
                else:
                    raise StandardError("At the moment, random_sample_fastq only supports paired end in this step config")
            

                # in-R1, in-R2, out-R1, out-R2
                if not len(input_paths) == 2:
                    raise StandardError("At the moment, exactly two input files are required.")

                input_files = misc.assign_strings(input_paths, ['R1', 'R2'])
                run.add_private_info('in-R1', input_files['R1'])
                run.add_private_info('in-R2', input_files['R2'])

    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')

            random_sample_fastq = [self.get_tool('random_sample_fastq'),
                                   self.get_option['readtype'], 
                                   '-z', 'y', '-n', str(self.get_option['reads']),
                                   '--FastqFileOut', 
                                   fifo_out_R1,
                                   fifo_out_R2,
                                   '--FastqFile',
                                   run.get_private_info('in-R1'),
                                   run.get_private_info('in-R2')]
                                   
            pool.launch(random_sample_fastq)

            pigz1 = [self.get_tool('pigz'), '--stdout', '--processes', '1',
                     '--blocksize', '4096', fifo_out_R1]
            pigz2 = [self.get_tool('pigz'), '--stdout', '--processes', '1',
                     '--blocksize', '4096', fifo_out_R2]
            output_files = run.get_output_files_for_annotation_and_tags('reads', ['R1', 'R2']))
            pool.launch(pigz1, stdout_path = run.get_private_info('out-R1'))
            pool.launch(pigz2, stdout_path = run.get_private_info('out-R2'))
            
