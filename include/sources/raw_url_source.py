import sys
from abstract_step import *
import os
import process_pool
import urlparse

class RawUrlSource(AbstractStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')
        self.require_tool('curl')
        self.require_tool('sha1sum')
        self.add_option('url', str, optional=False, description="file URL")
        self.add_option('sha1', str, optional=True, description="expected SHA1 checksum of downloaded file")
        
    def declare_runs(self):
        path = os.path.basename(urlparse.urlparse(self.option('url')).path)
        with self.declare_run('download') as run:
            run.add_output_file('raw', path, [])

    def execute(self, run_id, run):
        
        path = run.get_single_output_file_for_annotation('raw')
        
        # 1. download file and pipe to sha1sum
        download_sha1_path = self.get_temporary_path()
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                curl = [self.tool('curl'), self.option('url')]
                sha1sum = [self.tool('sha1sum'), '-b', '-']
                
                pipeline.append(curl, stdout_path = path)
                pipeline.append(sha1sum, stdout_path = download_sha1_path)
            
        if self.option_set_in_config('sha1'):
            with open(download_sha1_path, 'r') as f:
                line = f.read().strip().split(' ')
                if line[0] != self.option('sha1'):
                    # rename the output file, so the run won't be completed successfully
                    os.rename(path, path + '.mismatching.sha1')
                    raise StandardError("Error: SHA1 mismatch.")

        # remove the temporary file so that the temp directory can be deleted
        os.unlink(download_sha1_path)
