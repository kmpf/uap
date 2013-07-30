import sys
from abstract_step import *
import process_pool

class FixCutadapt(AbstractStep):
    
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)
        
        self.set_cores(6)

        self.add_connection('in/reads')
        self.add_connection('out/reads')

        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('fix_cutadapt')

    def declare_runs(self):
        # fetch all incoming run IDs which produce reads...
        run_ids = dict()
        for run_id, input_paths in self.run_ids_and_input_files_for_connection('in/reads'):
            fixed_run_id = run_id[:-3]
            which = run_id[-2:]
            if not fixed_run_id in run_ids:
                run_ids[fixed_run_id] = dict()
            if len(input_paths) != 1:
                raise StandardError("Expected one input file.")
            run_ids[fixed_run_id][which] = input_paths[0]
            
        for run_id, input_paths in run_ids.items():
            paired_end = self.find_upstream_info(run_id, 'paired_end')
            if not paired_end:
                raise StandardError("Not a paired end sample: %s" % run_id)
            
            with self.declare_run(run_id) as run:
                for which in ['R1', 'R2']:
                    out_path = "%s-fixed-%s.fastq.gz" % (run_id, which)
                    run.add_output_file("reads", out_path, input_paths.values())
                    run.add_private_info('in-%s' % which, input_paths[which])
            
    def execute(self, run_id, run):
        out_paths = run.get_output_files_for_annotation_and_tags('reads', ['R1', 'R2'])
        
        with process_pool.ProcessPool(self) as pool:
            fifo_in_R1 = pool.get_temporary_fifo('fifo_in_R1', 'input')
            fifo_in_R2 = pool.get_temporary_fifo('fifo_in_R2', 'input')
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), run.private_info('in-R1')]
                pigz = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = fifo_in_R1)
        
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), run.private_info('in-R2')]
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
                pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = out_paths['R1'])
                
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), fifo_out_R2]
                pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = out_paths['R2'])
