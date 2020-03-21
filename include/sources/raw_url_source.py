from uaperrors import StepError
import sys
from logging import getLogger
import os
import urllib.parse
from abstract_step import AbstractSourceStep

logger = getLogger("uap_logger")


class RawUrlSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')

        self.require_tool('compare_secure_hashes')
        # Step was tested for cp (GNU coreutils) release 8.25
        self.require_tool('cp')
        # Step was tested for curl release 7.47.0
        self.require_tool('curl')
        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        # Step was tested for mkdir (GNU coreutils) release 8.25
        self.require_tool('mkdir')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

        self.add_option('filename', str, optional=True,
                        description="local file name of downloaded file")
        self.add_option('hashing-algorithm', str, optional=True,
                        choices=['md5', 'sha1', 'sha224', 'sha256',
                                 'sha384', 'sha512'],
                        description="hashing algorithm to use")
        self.add_option('path', str, optional=False,
                        description="directory to move downloaded file to")
        self.add_option('secure-hash', str, optional=True,
                        description="expected secure hash of downloaded file")
        self.add_option('uncompress', bool, optional=True, default=False,
                        description='File is uncompressed after download')
        self.add_option('url', str, optional=False,
                        description="Download URL")

        # Options for dd
        self.add_option('dd-blocksize', str, optional=True, default="256k")

    def runs(self, run_ids_connections_files):
        # Get file name of downloaded file
        url_filename = os.path.basename(
            urllib.parse.urlparse(self.get_option('url')).path)

        # Is downloaded file gzipped?
        root, ext = os.path.splitext(url_filename)
        is_gzipped = True if ext in ['.gz', '.gzip'] else False
        if not is_gzipped and self.get_option('uncompress'):
            raise StepError(self, "Uncompression of non-gzipped file %s requested."
                           % url_filename)

        # Handle the filename to have the proper ending
        filename = root if self.get_option('uncompress') and is_gzipped \
            else url_filename

        if self.is_option_set_in_config('filename'):
            conf_filename = self.get_option('filename')
            root, ext = os.path.splitext(
                os.path.basename(conf_filename))

            if is_gzipped and self.get_option('uncompress') and \
               ext in ['.gz', '.gzip']:
                raise StepError(self, "The filename %s should NOT end on '.gz' or "
                               "'.gzip'." % conf_filename)
            filename = conf_filename

        # Get directory to move downloaded file to
        path = self.get_option('path')
        # Absolute path to downloaded file
        final_abspath = os.path.join(path, filename)

        with self.declare_run('download') as run:
            # Test if path exists
            if os.path.exists(path):
                # Fail if it is not a directory
                if not os.path.isdir(path):
                    raise StepError(self,
                        "Path %s already exists but is not a directory" % path)
            else:
                # Create the directory
                with run.new_exec_group() as mkdir_exec_group:
                    mkdir = [self.get_tool('mkdir'), '-p', path]
                    mkdir_exec_group.add_command(mkdir)
            out_file = run.add_output_file('raw', final_abspath, [])

            temp_filename = run.add_temporary_file(suffix=url_filename)
            with run.new_exec_group() as curl_exec_group:
                # 1. download file
                curl = [self.get_tool('curl'),
                        '--output', temp_filename,
                        self.get_option('url')]
                curl_exec_group.add_command(curl)

            if self.is_option_set_in_config('hashing-algorithm') and \
               self.is_option_set_in_config('secure-hash'):
                with run.new_exec_group() as check_exec_group:
                    # 2. Compare secure hashes
                    compare_secure_hashes = [
                        self.get_tool('compare_secure_hashes'),
                        '--algorithm',
                        self.get_option('hashing-algorithm'),
                        '--secure-hash',
                        self.get_option('secure-hash'),
                        temp_filename
                    ]
                    check_exec_group.add_command(compare_secure_hashes)
            with run.new_exec_group() as cp_exec_group:
                if self.get_option("uncompress"):
                    with cp_exec_group.add_pipeline() as pipe:
                        pigz = [self.get_tool('pigz'),
                                '--decompress',
                                '--stdout',
                                '--processes', '1',
                                temp_filename]
                        dd_out = [self.get_tool('dd'),
                                  'bs=%s' % self.get_option('dd-blocksize'),
                                  'of=%s' % out_file]
                        pipe.add_command(pigz)
                        pipe.add_command(dd_out)
                else:
                    cp = [self.get_tool('cp'), '--update', temp_filename,
                          out_file]
                    cp_exec_group.add_command(cp)
