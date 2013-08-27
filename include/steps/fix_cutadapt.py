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
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')

            if is_paired_end:
                # chop of '-R1' or '-R2' respectively
                fixed_run_id = run_id[:-3]
                # but keep info if 'R1' or 'R2'
                which = run_id[-2:]
                if not fixed_run_id in run_ids:
                    run_ids[fixed_run_id] = list()
                if len(input_paths) != 1:
                    raise StandardError("Expected one input file.")
                run_ids[fixed_run_id].extend(input_paths)
            else:
                if len(input_paths) != 1:
                    raise StandardError("Expected one input file.")
                run_ids[run_id] = input_paths

        for run_id, input_paths in run_ids.items():
            is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')

            if is_paired_end:
                with self.declare_run(run_id) as run:
                    run.add_private_info('paired_end', is_paired_end)
                    if len(input_paths) != 2:
                        raise StandardError("Expected two input files and got %s." % input_paths)
                    input_files = misc.assign_strings(input_paths, ['R1', 'R2'])
                    for which in ['R1', 'R2']:
                        out_path = "%s-fixed-%s.fastq.gz" % (run_id, which)
                        run.add_output_file("reads", out_path, input_paths)
                        run.add_private_info('in-%s' % which, input_files[which])          
            else:
                with self.declare_run(run_id) as run:
                    run.add_private_info('paired_end', is_paired_end)
                    out_path = "%s-fixed.fastq.gz" % run_id
                    run.add_output_file("reads", out_path, input_paths)
                    run.add_private_info('in-R1', input_paths)
            
    def execute(self, run_id, run):
        is_paired_end = run.get_private_info('paired_end')
        out_paths = None
        if is_paired_end:
            out_paths = run.get_output_files_for_annotation_and_tags('reads', ['R1', 'R2'])
        else:
            out_paths = run.get_output_files_for_annotation_and_tags('reads', [run_id])
        print(out_paths)
        with process_pool.ProcessPool(self) as pool:
            fifo_in_R1 = pool.get_temporary_fifo('fifo_in_R1', 'input')
            fifo_in_R2 = None
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            fifo_out_R2 = None
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
        
            fix_cutadapt = [self.get_tool('fix_cutadapt'), fifo_in_R1, fifo_out_R1]
            if is_paired_end:
                fix_cutadapt.extend(['--R2-in', fifo_in_R2])
            if is_paired_end:
                fix_cutadapt.extend(['--R2-out', fifo_out_R2])
            pool.launch(fix_cutadapt)
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), fifo_out_R1]
                pigz = [self.get_tool('pigz'), 
                        '--blocksize', '4096', 
                        '--processes', '2', 
                        '--stdout']
                pipeline.append(cat4m)
                out_path = None
                if is_paired_end:
                    out_path = out_paths['R1']
                else:
                    out_path = out_paths
                pipeline.append(pigz, stdout_path = out_path)
                
            if is_paired_end:
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), fifo_out_R2]
                    pigz = [self.get_tool('pigz'), 
                            '--blocksize', '4096', 
                            '--processes', '2', 
                            '--stdout']
                    pipeline.append(cat4m)
                    pipeline.append(pigz, stdout_path = out_paths['R2'])
