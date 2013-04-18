import sys
from abstract_step import *
import unix_pipeline
import tempfile
import yaml


class Segemehl(AbstractStep):
    def __init__(self, pipeline):
        super(Segemehl, self).__init__(pipeline)
        self.set_cores(12)

    def setup_runs(self, complete_input_run_info):
        # make sure tools are available
        self.tool('segemehl')
        self.tool('pigz')

        # make sure files are available
        for key in ['genome', 'index']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find " + key + " file: " + self.options[key])

        output_run_info = {}
        for run_id, input_run_info in complete_input_run_info.items():
            output_run_info[run_id] = {}
            output_run_info[run_id]['output_files'] = {}
            output_run_info[run_id]['output_files']['alignments']  = {}
            output_run_info[run_id]['output_files']['alignments'][run_id + '-segemehl-results.sam.gz'] = input_run_info['output_files']['reads'].keys()
            output_run_info[run_id]['output_files']['log']  = {}
            output_run_info[run_id]['output_files']['log'][run_id + '-segemehl-log.txt'] = input_run_info['output_files']['reads'].keys()
            read_files = assign_strings(input_run_info['output_files']['reads'].keys(), ['R1', 'R2'])
            output_run_info[run_id]['info'] = {}
            output_run_info[run_id]['info']['R1-in'] = read_files['R1']
            output_run_info[run_id]['info']['R2-in'] = read_files['R2']

        return output_run_info

    def execute(self, run_id, run_info):
        out_name = run_info['output_files']['alignments'].keys()[0]
        if out_name[-3:] != '.gz':
            raise StandardError("Expected .gz in output file name")

        temp_out_name = out_name[:-3]

        _, fifo_path_genome = tempfile.mkstemp('segemehl-genome-fifo')
        os.close(_)
        os.unlink(fifo_path_genome)
        os.mkfifo(fifo_path_genome)

        subprocess.Popen(
            [self.tool('cat4m'), self.options['genome'], '-o', fifo_path_genome],
            preexec_fn = os.setsid)

        segemehl = [
            self.tool('segemehl'),
            '-d', fifo_path_genome,
            '-i', self.options['index'],
            '-q', run_info['info']['R1-in'],
            '-p', run_info['info']['R2-in'],
            '-H', '1',
            '-t', '10',
            '-s', '-S',
            '-o', '/dev/stdout'
        ]
        
        pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '2', '-c']
        
        up = unix_pipeline.UnixPipeline()
        up.append(segemehl, stderr = open(run_info['output_files']['log'].keys()[0], 'w'))
        up.append(pigz, stdout = open(out_name, 'w'))
        up.run()

        os.unlink(fifo_path_genome)
