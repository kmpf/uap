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
        

        
        

    def setup_runs(self, complete_input_run_info, connection_info):
        
        
        print(yaml.dump(complete_input_run_info, default_flow_style = False))
        print(yaml.dump(connection_info, default_flow_style = False))



            
        output_run_info = {}
        '''
        RIB0001:
         ;kl d;lks g;kf g

         RIB00001: Run
        for run_id, step_info in connection_info['in/reads']['runs'].items():
            paired_end = self.find_upstream_info(run_id, 'paired_end')
            with self.create_run(run_id) as run:
            
                f = run.add_output_file(outpath)
                f.add_input_file(...)
                
                run.add_output_file(outpath, ( ... ))
                run.add_input_file_for_output_file(outpath, ...)
                
                run.set_public_info('paired_end', True)
                run.set_private_info('whoohoo', 'yesss')
        '''
                        
        for run_id, step_info in connection_info['in/reads']['runs'].items():
            paired_end = self.find_upstream_info(run_id, 'paired_end')

            if paired_end == True:
                out_path_R1 = "%s-subsample-R1.fastq.gz" % run_id
                out_path_R2 = "%s-subsample-R2.fastq.gz" % run_id
                output_run_info[run_id] = { 
                    'output_files': {
                        'reads': {
                            out_path_R1: [],
                            out_path_R2: []
                            }
                        },
                    'info': { },
                    'private_info': {
                        'paired_end': paired_end
                        }
                    }
            else:
                raise StandardError("At the moment, random_sample_fastq only supports paired end in this step config")
            
                # in-R1, in-R2, out-R1, out-R2
            
            if (len(step_info) != 1):
                raise StandardError("At the moment, random_sample_fastq may only have one parent which produces reads.")
            for step_name, paths in step_info.items():
                if len(paths) != 2:
                    raise StandardError("At the moment, exactly two input files are required.")
                x = misc.assign_strings(paths, ['R1', 'R2'])
                output_run_info[run_id]['info']['in-R1'] = x['R1']
                output_run_info[run_id]['info']['in-R2'] = x['R2']
                output_run_info[run_id]['info']['out-R1'] = out_path_R1
                output_run_info[run_id]['info']['out-R2'] = out_path_R2
                for path in paths:
                    output_run_info[run_id]['output_files']['reads'][out_path_R1].append(path)
                    output_run_info[run_id]['output_files']['reads'][out_path_R2].append(path)
                

        print(yaml.dump(output_run_info, default_flow_style = False))

        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')



            "chain for i in opt or sub fuction"

            random_sample_fastq = [self.get_tool('random_sample_fastq'),
                                   self.get_option['readtype'], 
                                   '-z', 'y', '-n', str(self.get_option['reads']),
                                   '--FastqFileOut', 
                                   fifo_out_R1,
                                   fifo_out_R2,
                                   '--FastqFile',
                                   run_info['info']['in-R1'],
                                   run_info['info']['in-R2']]
                                   
            pigz1 = [self.get_tool('pigz'), '--stdout', '--processes', '1', '--blocksize', '4096', fifo_out_R1]
            pigz2 = [self.get_tool('pigz'), '--stdout', '--processes', '1', '--blocksize', '4096', fifo_out_R2]

            pool.launch(random_sample_fastq)
            pool.launch(pigz1, stdout_path = run_info['info']['out-R1'])
            pool.launch(pigz2, stdout_path = run_info['info']['out-R2'])
            
