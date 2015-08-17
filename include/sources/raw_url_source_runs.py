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

        self.require_tool('compare_secure_hashes')
        self.require_tool('curl')
        self.require_tool('mv')

        self.add_option('filename', str, optional = True,
                        description = "local file name of downloaded file")
        self.add_option('hashing-algorithm', str, optional=False,
                        choices = ['md5', 'sha1', 'sha224', 'sha256',
                                   'sha384', 'sha512'],
                        description = "hashing algorithm to use")
        self.add_option('secure-hash', str, optional = False,
                        description = "expected secure hash of downloaded file")
        self.add_option('url', str, optional=False,
                        description = "file URL")
        
    def runs(self, run_ids_connections_files):
        # Get file name of downloaded file
        filename = os.path.basename(
            urlparse.urlparse(self.get_option('url')).path)
        if self.is_option_set_in_config('filename'):
            filename = self.get_option('filename')

        with self.declare_run('download') as run:
            out_file = run.add_output_file('raw', filename, [])
            temp_filename = run.add_temporary_file(filename)
            with run.new_exec_group() as curl_exec_group:
                # 1. download file
                curl = [self.get_tool('curl'), self.get_option('url')]
                curl_exec_group.add_command(curl, stdout_path = temp_filename)
            with run.new_exec_group() as check_exec_group:
                # 2. Compare secure hashes
                compare_secure_hashes = [self.get_tool('compare_secure_hashes'),
                                         '--algorithm',
                                         self.get_option('hashing-algorithm'),
                                         '--secure-hash',
                                         self.get_option('secure-hash'),
                                         temp_filename
                                     ]
                check_exec_group.add_command(compare_secure_hashes)
            with run.new_exec_group() as mv_exec_group:
                mv = [self.get_tool('mv'), temp_filename, out_file]
                mv_exec_group.add_command(mv)
