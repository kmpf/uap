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


    def setup_runs(self, complete_input_run_info, connection_info):
        # make sure files are available
        for key in ['genome', 'index']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find %s file: %s" % (key, self.options[key]))

        if not 'swap_reads' in self.options:
            self.options['swap_reads'] = False
            
        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for run_id, input_run_info in step_input_info.items():
                output_run_info[run_id] = {}
                output_run_info[run_id]['output_files'] = {}
                output_run_info[run_id]['output_files']['alignments']  = {}
                output_run_info[run_id]['output_files']['alignments'][run_id + '-segemehl-results.sam.gz'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['unmapped']  = {}
                output_run_info[run_id]['output_files']['unmapped'][run_id + '-segemehl-unmapped.fastq.gz'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['log']  = {}
                output_run_info[run_id]['output_files']['log'][run_id + '-segemehl-log.txt'] = input_run_info['output_files']['reads'].keys()
                read_files = misc.assign_strings(input_run_info['output_files']['reads'].keys(), ['R1', 'R2'])
                output_run_info[run_id]['info'] = {}
                output_run_info[run_id]['info']['R1-in'] = read_files['R1']
                output_run_info[run_id]['info']['R2-in'] = read_files['R2']

        return output_run_info

    def execute(self, run_id, run_info):
        
        out_name = run_info['output_files']['alignments'].keys()[0]
        if out_name[-3:] != '.gz':
            raise StandardError("Expected .gz in output file name")

        with process_pool.ProcessPool(self) as pool:
            
            fifo_path_genome = pool.get_temporary_fifo('segemehl-genome-fifo', 'input')
            fifo_path_unmapped = pool.get_temporary_fifo('segemehl-unmapped-fifo', 'output')
            
            pool.launch([self.tool('cat4m'), self.options['genome'], '-o', fifo_path_genome])
            
            with pool.Pipeline(pool) as pipeline:
            
                q = run_info['info']['R1-in']
                p = run_info['info']['R2-in']
                
                if self.options['swap_reads']:
                    q = run_info['info']['R2-in']
                    p = run_info['info']['R1-in']
                    
                segemehl = [
                    self.tool('segemehl'),
                    '-d', fifo_path_genome,
                    '-i', self.options['index'],
                    '-q', q,
                    '-p', p,
                    '-u', fifo_path_unmapped,
                    '-H', '1',
                    '-t', '5',
                    '-s', '-S',
                    '-D', '0',
                    '-o', '/dev/stdout'
                ]
                
                pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
                
                pipeline.append(segemehl, stderr_path = run_info['output_files']['log'].keys()[0])
                pipeline.append(pigz, stdout_path = out_name)
                
            with pool.Pipeline(pool) as pipeline:
                pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
                
                pipeline.append([self.tool('cat4m'), fifo_path_unmapped])
                pipeline.append(pigz, stdout_path = run_info['output_files']['unmapped'].keys()[0])
                
