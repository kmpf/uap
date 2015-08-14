import sys
from abstract_step import *
import os
import process_pool
import urlparse

logger = logging.getLogger("uap_logger")

class RawUrlSource(AbstractStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')

        self.require_tool('curl')
        self.require_tool('sha1sum')
        self.require_tool('md5sum')

        self.add_option('md5', str, optional=True,
                        description="expected MD5 checksum of downloaded file")
        self.add_option('filename', str, optional=True,
                        description="local file name of downloaded file")
        self.add_option('sha1', str, optional=True,
                        description="expected SHA1 checksum of downloaded file")
        self.add_option('url', str, optional=False,
                        description="file URL")
        
    def runs(self, run_ids_connections_files):
        # Get file name of downloaded file
        filename = os.path.basename(
            urlparse.urlparse(self.get_option('url')).path)
        if self.is_option_set_in_config('filename'):
            filename = self.get_option('filename')

        with self.declare_run('download') as run:
            run.add_output_file('raw', filename, [])
            with run.new_exec_group() as exec_group:
                with exec_group.add_pipeline() as curl_pipe:
                    # 1. download file and pipe to sha1sum
                    download_sha1_path = run.add_temporary_file()
                    curl = [self.get_tool('curl'), self.get_option('url')]
                    sha1sum = [self.get_tool('sha1sum'), '-b', '-']
                
                    curl_pipe.append(curl, stdout_path = filename)
                    curl_pipe.append(
                        sha1sum,
                        stdout_path = run.add_temporary_file())

                    # Separate script needs to be developed to do the SHA1/MD5
                    # checking!!!
            
        if self.is_option_set_in_config('sha1'):
            with open(download_sha1_path, 'r') as f:
                line = f.read().strip().split(' ')
                if line[0] != self.get_option('sha1'):
                    # rename the output file, so the run won't be completed successfully
                    os.rename(path, path + '.mismatching.sha1')
                    raise StandardError("Error: SHA1 mismatch.")

        # remove the temporary file so that the temp directory can be deleted
        os.unlink(download_sha1_path)
