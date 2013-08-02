import sys
from abstract_step import *
import process_pool
import yaml


class TestRealign(AbstractStep):

    def __init__(self, pipeline):
        super(TestRealign, self).__init__(pipeline)
        
        self.set_cores(8)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/splicesites')
        self.add_connection('out/transrealigned')
        self.add_connection('out/log')
        
        self.require_tool('testrealign')
        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cat4m')
        self.require_tool('grep')
        self.require_tool('invertGood')

    def setup_runs(self, complete_input_run_info, connection_info):
        # make sure files are available
        for key in ['genome']:
            if not os.path.exists(self.get_option(key)):
                raise StandardError("Could not find %s file: %s" % (key, self.get_option(key) ))

        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for run_id, input_run_info in step_input_info.items():
                output_run_info[run_id] = {}
                output_run_info[run_id]['output_files'] = {}
                # TODO: only depend on bam, not on bai
                # TODO: make this more generic (query language-like)
                output_run_info[run_id]['output_files']['alignments'] = {}
                output_run_info[run_id]['output_files']['alignments'][run_id + '-realigned.sam.gz'] = input_run_info['output_files']['alignments'].keys()
                output_run_info[run_id]['output_files']['splicesites']  = {}
                output_run_info[run_id]['output_files']['splicesites'][run_id + '-splicesites.bed'] = input_run_info['output_files']['alignments'].keys()
                output_run_info[run_id]['output_files']['transrealigned']  = {}
                output_run_info[run_id]['output_files']['transrealigned'][run_id + '-transrealigned.bed'] = input_run_info['output_files']['alignments'].keys()
                output_run_info[run_id]['output_files']['log']  = {}
                output_run_info[run_id]['output_files']['log'][run_id + '-testrealign-log.txt'] = input_run_info['output_files']['alignments'].keys()
                
                output_run_info[run_id]['info'] = {}
                output_run_info[run_id]['info']['bam-in'] = [_ for _ in input_run_info['output_files']['alignments'].keys() if _[-4:] == '.bam'][0]
                output_run_info[run_id]['info']['maxdist'] = '100'
                if self.get_option('maxdist'):
                    output_run_info[run_id]['info']['maxdist'] = str(self.get_option('maxdist') )


        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            
            fifo_path_genome = pool.get_temporary_fifo('genome-fifo', 'input')
            fifo_path_splicesites = pool.get_temporary_fifo('splicesites-fifo', 'output')
            fifo_path_transrealigned = pool.get_temporary_fifo('transrealigned-fifo', 'output')
            
            pool.launch([self.get_tool('cat4m'), self.get_option('genome'), '-o', fifo_path_genome])
            
            with pool.Pipeline(pool) as pipeline:
                cat4m = [
                    self.get_tool('cat4m'),
                    run_info['info']['bam-in']
                ]
                
                samtools = [
                    self.get_tool('samtools'),
                    'view', '-h', '-'
                ]
                
                grep = [
                    self.get_tool('grep'),
                    '-v',
                    "\t\\*\t"
                ]
                
                invertGood = [
                    self.get_tool('invertGood')
                ]

                testrealign = [
                    self.get_tool('testrealign'),
                    '-q', '/dev/stdin',
                    '-d', fifo_path_genome,
                    '-t', '4',
                    '-M', run_info['info']['maxdist'],
                    '-o', '/dev/stdout',
                    '-U', fifo_path_splicesites,
                    '-T', fifo_path_transrealigned
                ]
                
                pigz = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
                
                pipeline.append(cat4m)
                pipeline.append(samtools)
                pipeline.append(grep)
                pipeline.append(invertGood)
                pipeline.append(testrealign, stderr_path = run_info['output_files']['log'].keys()[0])
                pipeline.append(pigz, stdout_path = run_info['output_files']['alignments'].keys()[0])
            
            pool.launch([self.get_tool('cat4m'), fifo_path_splicesites], stdout_path = run_info['output_files']['splicesites'].keys()[0])
            
            pool.launch([self.get_tool('cat4m'), fifo_path_transrealigned], stdout_path = run_info['output_files']['transrealigned'].keys()[0])
