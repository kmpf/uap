import sys
from abstract_step import *
import process_pool

class RemoveReadsSegemehl(AbstractStep):
    
    def __init__(self, pipeline):
        super(RemoveReadsSegemehl, self).__init__(pipeline)
        
        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log')

        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('remove_reads_segemehl')
        self.require_tool('samtools')

    def setup_runs(self, complete_input_run_info, connection_info):
        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for input_run_id, input_run_info in step_input_info.items():
                new_key = input_run_id
                if not new_key in output_run_info:
                    output_run_info[new_key] = { 'output_files': { 'alignments': {}, 'log': {} }, 'info': {} }
                if len(input_run_info['output_files']['alignments'].keys()) != 1:
                    raise StandardError("Only one alignment input file allowed (fix this via constraints).")
                for in_path in sorted(input_run_info['output_files']['alignments'].keys()):
                    k2 = new_key + '-fixed.sam.gz'
                    if not k2 in output_run_info[new_key]['output_files']['alignments']:
                        output_run_info[new_key]['output_files']['alignments'][k2] = []
                    output_run_info[new_key]['output_files']['alignments'][k2].append(in_path)
                    k2_log = k2.replace('.sam.gz', '-log.txt')
                    if not k2_log in output_run_info[new_key]['output_files']['log']:
                        output_run_info[new_key]['output_files']['log'][k2_log] = []
                    output_run_info[new_key]['output_files']['log'][k2_log].append(in_path)
                    output_run_info[new_key]['info']['alignments-in'] = in_path
                    output_run_info[new_key]['info']['alignments-out'] = k2
                    output_run_info[new_key]['info']['log'] = k2_log
        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), run_info['info']['alignments-in']]
                pigz1 = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                samtools = [self.tool('samtools'), 'view', '-h', '-']                
                remove_reads_segemehl = [self.tool('remove_reads_segemehl')]
                pigz2 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2']
                
                pipeline.append(cat4m)
                if run_info['info']['alignments-in'][-7:] == '.sam.gz':
                    pipeline.append(pigz1)
                else:
                    pipeline.append(samtools)
                pipeline.append(remove_reads_segemehl, stderr_path = run_info['info']['log'])
                pipeline.append(pigz2, stdout_path = run_info['info']['alignments-out'])
