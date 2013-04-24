import sys
from abstract_step import *
import copy
import pipeline
import unix_pipeline
import yaml

class Head(AbstractStep):
    '''
    The head step filters the first few lines of any input file (1000 by 
    default). Uncompressed and Gzip-compressed files are handled correctly.
    
    Options:
    
    - ``lines``: specify the number of lines that should be returned
    '''
    
    def __init__(self, pipeline):
        super(Head, self).__init__(pipeline)
        self.set_cores(6)

    def setup_runs(self, input_run_info):

        # make sure tools are available
        self.tool('dd')
        self.tool('pigz')
        self.tool('head')

        count = 1000
        if 'lines' in self.options:
            count = self.options['lines']

        output_run_info = {}
        for input_run_id in input_run_info.keys():
            output_run_info[input_run_id] = {
                'output_files': {},
                'info': {}
            }
            if 'info' in input_run_info[input_run_id]:
                output_run_info[input_run_id]['info'] = input_run_info[input_run_id]['info']
            output_run_info[input_run_id]['info']['head-count'] = count
            for annotation, input_files in input_run_info[input_run_id]['output_files'].items():
                output_run_info[input_run_id]['output_files'][annotation] = {}
                for in_path in input_files.keys():
                    out_path = in_path.replace('.fastq.gz', '-head.fastq.gz')
                    out_path = copy.copy(in_path)
                    if '.' in out_path:
                        offset = out_path.index('.')
                        out_path = out_path[:offset] + '-head' + out_path[offset:]
                    else:
                        out_path = out_path + '-head'
                    output_run_info[input_run_id]['output_files'][annotation][out_path] = [in_path]
        return output_run_info

    def execute(self, run_id, run_info):
        for annotation in run_info['output_files'].keys():
            for outpath, inpaths in run_info['output_files'][annotation].items():
                if len(inpaths) != 1:
                    raise StandardError("Expected one input file per output file.")

                inpath = inpaths[0]
                count = run_info['info']['head-count']

                if inpath[-3:] == '.gz':
                    # set up processes for a gz-compressed file
                    pigz1 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '1',
                        '-d', '-c', inpath]

                    head = ['head', '-n', str(count)]

                    pigz2 = [self.tool('pigz'),
                        '--blocksize', '4096', '--processes', '3', '-c']

                    # create the pipeline and run it
                    up = unix_pipeline.create_pipeline()
                    up.append(pigz1)
                    up.append(head)
                    up.append(pigz2, stdout = open(outpath, 'w'))

                    unix_pipeline.wait()
                else:
                    # it's not a gz-compressed file
                    dd1 = [self.tool('dd'), 'bs=4M', 'if=' + inpath]

                    head = ['head', '-n', str(count)]

                    dd2 = [self.tool('dd'), 'bs=4M', 'of=' + outpath]

                    # create the pipeline and run it
                    up = unix_pipeline.create_pipeline()
                    up.append(dd1)
                    up.append(head)
                    up.append(dd2)

                    unix_pipeline.wait()