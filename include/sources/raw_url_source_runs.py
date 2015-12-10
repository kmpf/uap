import sys
from abstract_step import AbstractSourceStep
import logging
import os
import urlparse

logger = logging.getLogger("uap_logger")

class RawUrlSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')

        self.require_tool('compare_secure_hashes')
        self.require_tool('curl')
        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('pigz')

        self.add_option('filename', str, optional = True,
                        description = "local file name of downloaded file")
        self.add_option('hashing-algorithm', str, optional=False,
                        choices = ['md5', 'sha1', 'sha224', 'sha256',
                                   'sha384', 'sha512'],
                        description = "hashing algorithm to use")
        self.add_option('path', str, optional = False,
                        description = "directory to move downloaded file to")
        self.add_option('secure-hash', str, optional = False,
                        description = "expected secure hash of downloaded file")
        self.add_option('uncompress', bool, optional = True, default = False,
                        description = 'Shall the file be uncompressed after '
                        'downloading')
        self.add_option('url', str, optional = False,
                        description = "file URL")
        
    def runs(self, run_ids_connections_files):
        # Get file name of downloaded file
        url_filename = os.path.basename(
            urlparse.urlparse(self.get_option('url')).path)
        filename = url_filename
        if self.is_option_set_in_config('filename'):
            filename = self.get_option('filename')
        # Get directory to move downloaded file to
        path = self.get_option('path')
        # Absolute path to downloaded file
        final_abspath = os.path.join(path, filename)

        with self.declare_run('download') as run:
            # Test if path exists
            if os.path.exists(path):
                # Fail if it is not a directory
                if not os.path.isdir(path):
                    raise StandardError(
                        "Path %s already exists but is not a directory" % path)
            else:
                # Create the directory
                with run.new_exec_group() as mkdir_exec_group:
                    mkdir = [self.get_tool('mkdir'), '-p', path]
                    mkdir_exec_group.add_command(mkdir)
            out_file = run.add_output_file('raw', final_abspath, [] )
            temp_filename = run.add_temporary_file(suffix = url_filename)
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
                                         temp_filename]
                check_exec_group.add_command(compare_secure_hashes)
            with run.new_exec_group() as mv_exec_group:
                if self.get_option("uncompress"):
                    pigz = [self.get_tool('pigz'), '--decompress',
                            '--processes', '1', temp_filename]
                    temp_filename = os.path.splitext(temp_filename)[0]
                    mv_exec_group.add_command(pigz)

                mv = [self.get_tool('mv'), '--no-clobber', temp_filename,
                      out_file]
                mv_exec_group.add_command(mv)
