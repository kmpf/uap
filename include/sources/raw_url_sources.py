from uaperrors import UAPError
import sys
from logging import getLogger
import os
import urllib.parse
from abstract_step import AbstractSourceStep

logger=getLogger("uap_logger")

class RawUrlSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(RawUrlSource, self).__init__(pipeline)

        self.add_connection('out/raw')

        self.require_tool('compare_secure_hashes')
        self.require_tool('cp')
        self.require_tool('curl')
        self.require_tool('dd')
        self.require_tool('pigz')

        self.add_option('run-download-info', dict, optional = False,
                        description = "Dictionary of dictionaries. "
                        "The keys are the names of the runs. "
                        "The values are dictionaries whose keys are identical "
                        "with the options of an 'raw_url_source' source step. "
                        "An example: "
                        "<name>:\n"
                        "    filename: <filename>\n"
                        "    hashing-algorithm: <hashing-algorithm>\n"
                        "    path: <path>\n"
                        "    secure-hash: <secure-hash>\n"
                        "    uncompress: <uncompress>\n"
                        "    url: <url>")

        # Options for dd
        self.add_option('dd-blocksize', str, optional = True, default = "256k")

    def runs(self, run_ids_connections_files):
        # Sanity check the 'file-download-map'
        file_download = self.get_option('run-download-info')

        # List with valid download options
        download_opts = {'filename', 'hashing-algorithm', 'path',
                            'secure-hash', 'uncompress', 'url'}
        mandatory_opts = {'filename', 'path', 'url'}

        for files, downloads in file_download.items():
            # Control input for unknown options
            unknown_opts = set( downloads.keys() ).difference(download_opts)
            if len( unknown_opts ) > 0:
                raise UAPError("Unknown option(s) %s for download of %s"
                             % (" ".join(unknown_opts), files) )
            # Control input for missing mandatory options
            missing_mandatory_opts = mandatory_opts.difference( set(
                downloads.keys() ))
            if len(missing_mandatory_opts) > 0:
                raise UAPError("Download of %s misses mandatory option(s): %s"
                             % (files, " ".join(missing_mandatory_opts)) )

            # Check the optional parameters and set default if not available
            downloads['hashing-algorithm'] = downloads.get('hashing-algorithm')
            downloads['secure-hash'] = downloads.get('secure-hash')
            downloads['uncompress'] = downloads.get('uncompress', False)

            # 1. Check the 'hashing-algorithm'
            hash_algos = ['md5', 'sha1', 'sha224', 'sha256', 'sha384', 'sha512']
            if downloads['hashing-algorithm'] not in hash_algos:
                raise UAPError("Option 'hashing-algorithm' for download %s "
                             "has invalid value %s. Has to be one of %s."
                             % (files, downloads['hashing-algorithm'],
                                ", ".join(hash_algos)) )

            # 2. Check the 'secure-hash'
            if isinstance(downloads['secure-hash'], str) and not \
                   downloads['hashing-algorithm']:
                    raise UAPError("Option 'secure-hash' set for download %s "
                                 "but option 'hashing-algorithm' is missing."
                                 % files)

            # Get file name of downloaded file
            url_filename = os.path.basename(
                urllib.parse.urlparse(downloads['url']).path)

            # Is downloaded file gzipped?
            root, ext = os.path.splitext(url_filename)
            is_gzipped = True if ext in ['.gz', '.gzip'] else False
            if not is_gzipped and downloads['uncompress']:
                raise UAPError("Uncompression of non-gzipped file %s requested."
                                    % url_filename)
            # Handle the filename to have the proper ending
            filename = root if downloads['uncompress'] and is_gzipped \
                       else url_filename

            conf_filename = downloads['filename']
            root, ext = os.path.splitext(
                os.path.basename(conf_filename))

            if is_gzipped and downloads['uncompress'] and \
               ext in ['.gz', '.gzip']:
                raise UAPError("The filename %s should NOT end on '.gz' or "
                             "'.gzip'." % conf_filename)
            filename = conf_filename

            # Get directory to move downloaded file to
            path = downloads['path']
            # Absolute path to downloaded file
            final_abspath = os.path.join(path, filename)

            with self.declare_run(files) as run:
                # Test if path exists
                if os.path.exists(path):
                    # Fail if it is not a directory
                    if not os.path.isdir(path):
                        raise UAPError(
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
                    curl = [self.get_tool('curl'), downloads['url']]
                    curl_exec_group.add_command(curl, stdout_path = temp_filename)
            
                if downloads['hashing-algorithm'] and downloads['secure-hash']:
                    with run.new_exec_group() as check_exec_group:
                        # 2. Compare secure hashes
                        compare_secure_hashes = [
                            self.get_tool('compare_secure_hashes'),
                            '--algorithm',
                            downloads['hashing-algorithm'],
                            '--secure-hash',
                            downloads['secure-hash'],
                            temp_filename
                        ]
                        check_exec_group.add_command(compare_secure_hashes)

                with run.new_exec_group() as cp_exec_group:
                    if downloads["uncompress"]:
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
