import sys
from abstract_step import *
import misc
import unix_pipeline
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


    def setup_runs(self, complete_input_run_info):
        # make sure files are available
        for key in ['genome', 'index']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find %s file: %s" % (key, self.options[key]))

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

        fifo_path_genome = unix_pipeline.mkfifo('segemehl-genome-fifo')
        fifo_path_unmapped = unix_pipeline.mkfifo('segemehl-unmapped-fifo')
        
        unix_pipeline.launch([self.tool('cat4m'), self.options['genome'], '-o', fifo_path_genome])

        segemehl = [
            self.tool('segemehl'),
            '-d', fifo_path_genome,
            '-i', self.options['index'],
            '-q', run_info['info']['R1-in'],
            '-p', run_info['info']['R2-in'],
            '-u', fifo_path_unmapped,
            '-H', '1',
            '-t', '10',
            '-s', '-S',
            '-D', '0',
            '-o', '/dev/stdout'
        ]
        
        pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
        
        p = unix_pipeline.create_pipeline()
        p.append(segemehl, stderr = open(run_info['output_files']['log'].keys()[0], 'w'))
        p.append(pigz, stdout = open(out_name, 'w'))
        
        p2 = unix_pipeline.create_pipeline()
        pigz2 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
        p2.append([self.tool('cat4m'), fifo_path_unmapped])
        p2.append(pigz2, stdout = open(run_info['output_files']['unmapped'].keys()[0], 'w'))
        
        unix_pipeline.wait()

        os.unlink(fifo_path_genome)
        os.unlink(fifo_path_unmapped)
