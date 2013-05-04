import sys
from abstract_step import *
import unix_pipeline

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
        fifo_in_R1 = self.get_temporary_fifo('fifo_in_R1', 'input')
        fifo_in_R2 = self.get_temporary_fifo('fifo_in_R2', 'input')
        fifo_out_R1 = self.get_temporary_fifo('fifo_out_R1', 'output')
        fifo_out_R2 = self.get_temporary_fifo('fifo_out_R2', 'output')
        
        cat4m1 = [self.tool('cat4m'), run_info['info']['R1-in']]
        pigz1 = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
        
        p1 = unix_pipeline.UnixPipeline()
        p1.append(cat4m1)
        p1.append(pigz1, stdout_path = fifo_in_R1)
        
        cat4m2 = [self.tool('cat4m'), run_info['info']['R2-in']]
        pigz2 = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
        
        p2 = unix_pipeline.UnixPipeline()
        p2.append(cat4m2)
        p2.append(pigz2, stdout_path = fifo_in_R2)

        fix_cutadapt = [
            self.tool('fix_cutadapt'),
            fifo_in_R1,
            fifo_in_R2,
            fifo_out_R1,
            fifo_out_R2
        ]
        
        unix_pipeline.launch(fix_cutadapt)
        
        p3 = unix_pipeline.UnixPipeline()
        cat4m3 = [self.tool('cat4m'), fifo_out_R1]
        pigz3 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '1', '--stdout']
        
        p4 = unix_pipeline.UnixPipeline()
        cat4m4 = [self.tool('cat4m'), fifo_out_R2]
        pigz4 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '1', '--stdout']

        p3.append(cat4m3)
        p3.append(pigz3, stdout_path = run_info['info']['R1-out'])
        
        p4.append(cat4m4)
        p4.append(pigz4, stdout_path = run_info['info']['R2-out'])
        
        unix_pipeline.wait()
        
        os.unlink(fifo_in_R1)
        os.unlink(fifo_in_R2)
        os.unlink(fifo_out_R1)
        os.unlink(fifo_out_R2)
        