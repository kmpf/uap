import sys
from abstract_step import *
import os
import urlparse
import process_pool

class RawUrlSource(AbstractStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')
        self.require_tool('curl')
        self.add_option('url', str, optional=False, description="download url")
        self.add_option('sha1', str, optional=True, description="sha1 cheksum of file")
        
    def setup_runs(self, input_run_info, connection_info):
        output_run_info = {}
       
        path = os.path.basename(urlparse.urlparse(self.option('url')).path)
        output_run_info['download'] = {}
        output_run_info['download']['output_files'] = {}
        output_run_info['download']['output_files']['raw'] = {}
        output_run_info['download']['output_files']['raw'][path] = []

        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            curl = [self.tool('curl'), self.option('url')]
            pool.launch(curl, stdout_path = run_info['output_files']['raw'].keys()[0])
            
        # TODO: verify checksum after process pool has finished
        #stdout_assert_sha1 = self.options['sha1']
