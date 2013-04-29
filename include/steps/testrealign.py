import sys
from abstract_step import *
import unix_pipeline
import yaml


class TestRealign(AbstractStep):

    def __init__(self, pipeline):
        super(TestRealign, self).__init__(pipeline)
        
        self.set_cores(12)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/splicesites')
        self.add_connection('out/transrealigned')
        self.add_connection('out/log')
        
        self.require_tool('testrealign')
        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cat4m')

    def setup_runs(self, complete_input_run_info):
        # make sure files are available
        for key in ['genome']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find %s file: %s" % (key, self.options[key]))

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
                if 'maxdist' in self.options:
                    output_run_info[run_id]['info']['maxdist'] = str(self.options['maxdist'])


        return output_run_info

    def execute(self, run_id, run_info):
        fifo_path_genome = unix_pipeline.mkfifo('realigner-genome-fifo')
        fifo_path_splicesites = unix_pipeline.mkfifo('realigner-splicesites-fifo')
        fifo_path_transrealigned = unix_pipeline.mkfifo('realigner-transrealigned-fifo')
        
        unix_pipeline.launch([self.tool('cat4m'), self.options['genome'], '-o', fifo_path_genome])
        
        cat4m = [
            self.tool('cat4m'),
            run_info['info']['bam-in']
        ]
        
        samtools = [
            self.tool('samtools'),
            'view', '-h', '-'
        ]

        testrealign = [
            self.tool('testrealign'),
            '-q', '/dev/stdin',
            '-d', fifo_path_genome,
            '-t', '10',
            '-M', run_info['info']['maxdist'],
            '-o', '/dev/stdout',
            '-U', fifo_path_splicesites,
            '-T', fifo_path_transrealigned
        ]
        
        pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
        
        p = unix_pipeline.UnixPipeline()
        p.append(cat4m)
        p.append(samtools)
        p.append(testrealign, stderr_path = run_info['output_files']['log'].keys()[0])
        p.append(pigz, stdout_path = run_info['output_files']['alignments'].keys()[0])
        
        unix_pipeline.launch([self.tool('cat4m'), fifo_path_splicesites, '-o', run_info['output_files']['splicesites'].keys()[0]])
        unix_pipeline.launch([self.tool('cat4m'), fifo_path_transrealigned, '-o', run_info['output_files']['transrealigned'].keys()[0]])
        
        unix_pipeline.wait()

        os.unlink(fifo_path_genome)
        os.unlink(fifo_path_splicesites)
