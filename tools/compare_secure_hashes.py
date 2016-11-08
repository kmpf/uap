#!/usr/bin/env python
import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import argparse
import sys
import hashlib

def main():

    parser = argparse.ArgumentParser(
        description = 'This script computes a secure hash and compares it ' +
        'with a given hash value. If they do not fit the input file will be ' +
        'renamed and a error message is output.',
        prog = 'compare_checksums.py',
        formatter_class = argparse.RawTextHelpFormatter)

    parser.add_argument('--version',
                        action = 'version',
                        version = '%(prog)s 0.01'
    )

    parser.add_argument("file_to_hash",
                        help = "file to compare checksums for",
                        type = argparse.FileType('r')
    )

    parser.add_argument("--algorithm",
                        choices = ['md5', 'sha1', 'sha224', 'sha256',
                                   'sha384', 'sha512'],
                        default = None,
                        dest = "hash_alg",
                        required = True,
                        type = str,
                        help = "hashing algorithm to use"
    )

    parser.add_argument("--secure-hash",
                        dest = "provided_hash_value",
                        required = True,
                        type = str,
                        help = "secure hash used for comparision")

    # get arguments and call the appropriate function
    args = parser.parse_args()


    computed_hash_value = hashfile(args.file_to_hash,
                                   getattr(hashlib, args.hash_alg)()
                               )
    print("Provided hash value: %s" % args.provided_hash_value)
    print("Computed hash value: %s" % computed_hash_value)

    file_to_hash_abspath = os.path.abspath(args.file_to_hash.name)
    abspath, filename = os.path.split(file_to_hash_abspath)
    file_to_hash_new_path = os.path.join(
        abspath, "%s.mismatching.%s" % (filename, args.hash_alg))
    if args.provided_hash_value == computed_hash_value:
        sys.exit()
    else:
        if os.path.exists(file_to_hash_new_path):
            raise StandardError("File %s already exists. Couldn't rename %s "
                                "to %s." % (file_to_hash_new_path,
                                            file_to_hash_abspath,
                                            file_to_hash_new_path)
            )
        os.rename(file_to_hash_abspath, file_to_hash_new_path)
        sys.exit("Mismatching secure hashes! File %s was renamed to %s" %
                 (file_to_hash_abspath, file_to_hash_new_path) )

# Copied from http://stackoverflow.com/questions/3431825/generating-a-md5-checksum-of-a-file
def hashfile(file_handle, hasher, blocksize=65536):
    buf = file_handle.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = file_handle.read(blocksize)
    return hasher.hexdigest()

if __name__ == '__main__':
    main()
