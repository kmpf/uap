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
    This step is useful to drastically reduce the amount of data in order
    to quickly test a pipeline.
    
    Options:
    
    - ``lines``: specify the number of lines that should be returned
    '''
    
    def __init__(self, pipeline):
        super(Head, self).__init__(pipeline)
        
        self.set_cores(6)
        
        self.add_connection('out/*')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('head')

    def setup_runs(self, input_run_info, connection_info):
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
            for tag, input_files in input_run_info[input_run_id]['output_files'].items():
                output_run_info[input_run_id]['output_files'][tag] = {}
                for in_path in input_files.keys():
                    out_path = in_path.replace('.fastq.gz', '-head.fastq.gz')
                    out_path = copy.copy(in_path)
                    if '.' in out_path:
                        offset = out_path.index('.')
                        out_path = out_path[:offset] + '-head' + out_path[offset:]
                    else:
                        out_path = out_path + '-head'
                    output_run_info[input_run_id]['output_files'][tag][out_path] = [in_path]
        return output_run_info

    def execute(self, run_id, run_info):
        for tag in run_info['output_files'].keys():
            for outpath, inpaths in run_info['output_files'][tag].items():
                if len(inpaths) != 1:
                    raise StandardError("Expected one input file per output file.")

                inpath = inpaths[0]
                count = run_info['info']['head-count']

                if inpath[-3:] == '.gz':
                    # set up processes for a gz-compressed file
                    cat4m = [self.tool('cat4m'), inpath]
                    pigz1 = [self.tool('pigz'), '--processes', '2', '--decompress', '--stdout', '-']
                    head = ['head', '-n', str(count)]
                    pigz2 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '3', '-c']

                    # create the pipeline and run it
                    up = unix_pipeline.UnixPipeline()
                    up.append(cat4m)
                    up.append(pigz1)
                    up.append(head)
                    up.append(pigz2, stdout_path = outpath)

                    unix_pipeline.wait()
                else:
                    # it's not a gz-compressed file
                    cat4m = [self.tool('cat4m'), inpath]

                    head = ['head', '-n', str(count)]

                    # create the pipeline and run it
                    up = unix_pipeline.UnixPipeline()
                    up.append(cat4m)
                    up.append(head, stdout_path = outpath)

                    unix_pipeline.wait()
