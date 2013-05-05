import sys
from abstract_step import *
import process_pool

class FixCutadapt(AbstractStep):
    
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)
        
        self.set_cores(9)

        self.add_connection('in/reads')
        self.add_connection('out/reads')

        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('fix_cutadapt')

    def setup_runs(self, complete_input_run_info, connection_info):
        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for input_run_id, input_run_info in step_input_info.items():
                if not 'read_number' in input_run_info['info']:
                    raise StandardError("fix_cutadapt can only be run on paired-end sequenced samples.")
                new_key = input_run_id.replace('-R1', '').replace('-R2', '')
                if not new_key in output_run_info:
                    output_run_info[new_key] = { 'output_files': { 'reads': {} }, 'info': {} }
                output_run_info[new_key]['info'][input_run_info['info']['read_number'] + '-in'] = input_run_info['output_files']['reads'].keys()[0]
                for in_path in sorted(input_run_info['output_files']['reads'].keys()):
                    for _ in ['R1', 'R2']:
                        k2 = new_key + '-fixed-' + _ + '.fastq.gz'
                        if not k2 in output_run_info[new_key]['output_files']['reads']:
                            output_run_info[new_key]['output_files']['reads'][k2] = []
                        output_run_info[new_key]['output_files']['reads'][k2].append(in_path)
                        output_run_info[new_key]['info'][_ + '-out'] = k2
        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            fifo_in_R1 = pool.get_temporary_fifo('fifo_in_R1', 'input')
            fifo_in_R2 = pool.get_temporary_fifo('fifo_in_R2', 'input')
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), run_info['info']['R1-in']]
                pigz = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = fifo_in_R1)
        
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), run_info['info']['R2-in']]
                pigz = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = fifo_in_R2)
        
            fix_cutadapt = [
                self.tool('fix_cutadapt'),
                fifo_in_R1,
                fifo_in_R2,
                fifo_out_R1,
                fifo_out_R2
            ]
            
            pool.launch(fix_cutadapt)
        
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), fifo_out_R1]
                pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '1', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = run_info['info']['R1-out'])
                
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), fifo_out_R2]
                pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '1', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = run_info['info']['R2-out'])
