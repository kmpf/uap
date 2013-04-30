import sys
from abstract_step import *
import os
import urlparse
import unix_pipeline

class RawUrlSource(AbstractStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')
        self.require_tool('curl')

    def setup_runs(self, input_run_info, connection_info):
        output_run_info = {}

        if not 'url' in self.options:
            raise StandardError("missing 'url' key in raw_url_source")
        if not 'sha1' in self.options:
            raise StandardError("missing 'sha1' key in raw_url_source")
        
        
        path = os.path.basename(urlparse.urlparse(self.options['url']).path)
        output_run_info['download'] = {}
        output_run_info['download']['output_files'] = {}
        output_run_info['download']['output_files']['raw'] = {}
        output_run_info['download']['output_files']['raw'][path] = []

        return output_run_info

    def execute(self, run_id, run_info):
        # set up processes
        curl = [self.tool('curl'), self.options['url']]

        # create the pipeline and run it
        p = unix_pipeline.UnixPipeline()
        p.append(curl, stdout_path = run_info['output_files']['raw'].keys()[0])

        unix_pipeline.wait()
        
        if unix_pipeline.sha1_checksum_for_file_basename[os.path.basename(run_info['output_files']['raw'].keys()[0])] != self.options['sha1']:
            raise StandardError("Checksum of downloaded file does not match.")
