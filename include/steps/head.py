import sys
from abstract_step import *
import pipeline
import unix_pipeline
import yaml

class Head(AbstractStep):
    def __init__(self, pipeline):
        super(Head, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        #print(yaml.dump(input_run_info, default_flow_style = False))
        output_run_info = {}
        for input_run_id in input_run_info.keys():
            output_run_info[input_run_id] = { 'output_files': {} }
            for annotation, input_files in input_run_info[input_run_id]['output_files'].items():
                output_run_info[input_run_id]['output_files'][annotation] = {}
                for in_path in input_files.keys():
                    out_path = in_path.replace('.fastq.gz', '-head.fastq.gz')
                    output_run_info[input_run_id]['output_files'][annotation][out_path] = [in_path]
        return output_run_info

    def execute(self, run_id, run_info):
        for annotation in run_info['output_files'].keys():
            for outpath, inpaths in run_info['output_files'][annotation].items():
                if len(inpaths) != 1:
                    raise StandardError("Expected one input file per output file.")
                inpath = inpaths[0]

                # set up processes
                pigz1 = [self.pipeline.config['tools']['pigz']['path'], '-d', '-c', inpath]

                head = ['head', '-n', '1000']

                pigz2 = [self.pipeline.config['tools']['pigz']['path'],
                    '--blocksize', '4096', '--processes', '3', '-c']

                # create the pipeline and run it
                up = unix_pipeline.UnixPipeline()
                up.append(pigz1)
                up.append(head)
                up.append(pigz2, stdout = open(outpath, 'w'))

                up.run()
