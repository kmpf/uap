import sys
from abstract_step import *
import process_pool

class FixCutadapt(AbstractStep):
    
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)
        
        self.set_cores(6)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')

        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('fix_cutadapt')

    def declare_runs(self):
        # fetch all incoming run IDs which produce reads...
        found_files = dict()
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        paired_end_info = dict()

        for read in read_types.keys():
            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/%s' % read):
                if input_paths != [None]:
                    paired_end_info[run_id] = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')

                # save information per file in found_files
                if not run_id in found_files:
                    found_files[run_id] = dict()

                if not read in found_files[run_id]:
                    found_files[run_id][read] = list()
                # Check if we get exactly one input file
                if len(input_paths) != 1:
                    raise StandardError("Expected one input file.")
                found_files[run_id][read].extend(input_paths)

        for run_id in found_files.keys():
            with self.declare_run(run_id) as run:
                run.add_private_info('paired_end', paired_end_info[run_id])
                for read in found_files[run_id].keys():
                    run.add_output_file(read, 
                        "%s-fixed%s.fastq.gz" % (run_id, read_types[read]), 
                        found_files[run_id][read])
            
    def execute(self, run_id, run):
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        is_paired_end = run.get_private_info('paired_end')

        with process_pool.ProcessPool(self) as pool:
            first_read_path = run.get_single_output_file_for_annotation(
                'first_read')
            second_read_path = None
            fifo_in_R1 = pool.get_temporary_fifo('fifo_in_R1', 'input')
            fifo_in_R2 = None
            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
            fifo_out_R2 = None
            if is_paired_end:
                second_read_path = run.get_single_output_file_for_annotation(
                    'second_read')
                fifo_in_R2 = pool.get_temporary_fifo('fifo_in_R2', 'input')
                fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), run.get_input_files_for_output_file(first_read_path)]
                pigz = [self.get_tool('pigz'), '--decompress', '--processes', 
                        '1', '--stdout']                
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = fifo_in_R1)
                
            if is_paired_end:    
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), 
                             run.get_input_files_for_output_file(
                                 second_read_path)]
                    pigz = [self.get_tool('pigz'), '--decompress', 
                            '--processes', '1', '--stdout']
                    print(cat4m)
                    print(pigz)
                    print(fifo_in_R2)
                    pipeline.append(cat4m)
                    pipeline.append(pigz, stdout_path = fifo_in_R2)
        
            fix_cutadapt = [self.get_tool('fix_cutadapt'), fifo_in_R1, 
                            fifo_out_R1]
            if is_paired_end:
                fix_cutadapt.extend(['--R2-in', fifo_in_R2,
                                     '--R2-out', fifo_out_R2])
            pool.launch(fix_cutadapt)
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.get_tool('cat4m'), fifo_out_R1]
                pigz = [self.get_tool('pigz'), 
                        '--blocksize', '4096', 
                        '--processes', '2', 
                        '--stdout']
                pipeline.append(cat4m)
                pipeline.append(pigz, stdout_path = first_read_path)
                
            if is_paired_end:
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), fifo_out_R2]
                    pigz = [self.get_tool('pigz'), 
                            '--blocksize', '4096', 
                            '--processes', '2', 
                            '--stdout']
                    pipeline.append(cat4m)
                    pipeline.append(pigz, stdout_path = second_read_path)
