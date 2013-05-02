import sys
from abstract_step import *
import pipeline
import re
import unix_pipeline
import yaml

class HtSeqCount(AbstractStep):
    
    def __init__(self, pipeline):
        super(HtSeqCount, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/raw', constraints = {'total_files': 1})
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('grep')

    def setup_runs(self, complete_input_run_info, connection_info):
        
        output_run_info = {}
        
        if not 'pattern' in self.options:
            raise StandardError("pattern missing")
        
        for step_id, step_info in complete_input_run_info.items():
            for run_id, run_info in step_info.items():
                output_run_info[run_id] = dict()
                output_run_info[run_id]['info'] = dict()
                output_run_info[run_id]['output_files'] = dict()
                for tag, output_files in run_info['output_files'].items():
                    self.add_connection('out/' + tag)
                    output_run_info[run_id]['output_files'][tag] = dict()
                    for output_path, input_paths in output_files.items():
                        output_run_info[run_id]['info']['out_path'] = 'grepped-' + output_path
                        output_run_info[run_id]['info']['in_path'] = output_path
                        output_run_info[run_id]['output_files'][tag]['grepped-' + output_path] = [output_path]
        
        return output_run_info
    
    
    def execute(self, run_id, run_info):
        p = unix_pipeline.UnixPipeline()
        
        pigz1 = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout', run_info['info']['in_path']]
        grep = [self.tool('grep'), self.options['pattern']]
        pigz2 = [self.tool('pigz'), '--processes', '2', '--blocksize', '4096', '--stdout']
        
        p.append(pigz1)
        p.append(grep)
        p.append(pigz2, stdout_path = run_info['info']['out_path'])
                
        unix_pipeline.wait()
        