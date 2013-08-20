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
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
            print(type(input_paths))
            print(input_paths)
            if is_paired_end:
                with self.declare_run(run_id) as run:
                    run.add_private_info('paired_end', is_paired_end)
                    in_files = {'R1': list(), 'R2': list()}
                    for input_path in input_paths:
                        which = misc.assign_string(input_path, ['R1', 'R2'])
                        in_files[which].append(input_path)

                    for which in ['R1', 'R2']:
                        out_path = "%s-fixed-%s.fastq.gz" % (run_id, which)
                        run.add_output_file("reads", out_path, in_files.values())
                        run.add_private_info('in-%s' % which, in_files[which])
 
            else:
                with self.declare_run(run_id) as run:
                    run.add_private_info('paired_end', is_paired_end)
                    out_path = "%s-fixed.fastq.gz" % run_id
                    run.add_output_file("reads", out_path, input_paths)
                    run.add_private_info('in-R1', input_paths)
            
    def execute(self, run_id, run):
        is_paired_end = run.get_private_info('paired_end')
        if is_paired_end:
            out_paths = run.get_output_files_for_annotation_and_tags('reads', ['R1', 'R2'])
        else:
            out_paths = run.get_output_files_for_annotation_and_tags('reads', [run_id])
        with process_pool.ProcessPool(self) as pool:
            fifo_in_R1 = pool.get_temporary_fifo('fifo_in_R1', 'input')
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            if is_paired_end:
                fifo_in_R2 = pool.get_temporary_fifo('fifo_in_R2', 'input')
                fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), run.get_private_info('in-R1')]
                pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = fifo_in_R1)
                
            if is_paired_end:    
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), run.get_private_info('in-R2')]
                    pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                
                    pipeline.append(cat4m)
                    pipeline.append(pigz, stdout_path = fifo_in_R2)
        
            fix_cutadapt = [
                self.get_tool('fix_cutadapt'),
                fifo_in_R1
                ]
            if is_paired_end:
                fix_cutadapt.append(fifo_in_R2)
            fix_cutadapt.append(fifo_out_R1)
            if is_paired_end:
                fix_cutadapt.append(fifo_out_R2)
            
            pool.launch(fix_cutadapt)
        
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), fifo_out_R1]
                pigz = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '2', '--stdout']
                
                pipeline.append(cat4m)
                out_path = str()
                if is_paired_end:
                    out_path = out_paths['R1']
                else:
                    out_path = out_paths
                pipeline.append(pigz, stdout_path = out_path)
                
            if is_paired_end:
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), fifo_out_R2]
                    pigz = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '2', '--stdout']
                
                    pipeline.append(cat4m)
                    pipeline.append(pigz, stdout_path = out_paths['R2'])
