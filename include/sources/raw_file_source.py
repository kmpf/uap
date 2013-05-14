import sys
from abstract_step import *
import os
import urlparse
import process_pool

class RawFileSource(AbstractStep):

    def __init__(self, pipeline):
        super(RawFileSource, self).__init__(pipeline)

        self.add_connection('out/raw')
        self.require_tool('cat4m')

    def setup_runs(self, input_run_info, connection_info):
        output_run_info = {}

        if not 'path' in self.options:
            raise StandardError("missing 'path' key in raw_file_source")
        if not 'sha1' in self.options:
            raise StandardError("missing 'sha1' key in raw_file_source")
        
        path = os.path.basename(urlparse.urlparse(self.options['path']).path)
        output_run_info[''] = {}
        output_run_info['download']['output_files'] = {}
        output_run_info['download']['output_files']['raw'] = {}
        output_run_info['download']['output_files']['raw'][path] = []

        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            cat4m = [self.tool('cat4m'), self.options['path']]
            pool.launch(curl, stdout_path = run_info['output_files']['raw'].keys()[0])
            
        # TODO: verify checksum after process pool has finished
        #stdout_assert_sha1 = self.options['sha1']