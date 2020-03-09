#!/bin/bash
import sys
import argparse
"exec" "`dirname $0`/../python_env/bin/python" "$0" "$@"

# ^^^
# the cmd above ensures that the correct python environment is
# selected to execute this script.
# The correct environment is the one belonging to uap, since all
# neccessary python modules are installed there.


import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
exec(
    compile(
        open(activate_this_file).read(),
        activate_this_file,
        'exec'),
    dict(
        __file__=activate_this_file))


def main():

    # Definition of the argument parser

    parser = argparse.ArgumentParser(
        description="This script is able to repair non-standard QNAMEs in FASTQ"
        " or SAM files.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "input_file",
        nargs='?',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help="Name of input file where QNAMES need to be "
        "fixed. Reads by default from STDIN.")

    parser.add_argument(
        "output_file",
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Name of output file to which corrected data is written to. "
        "Reads by default from STDOUT.",
    )

    parser.add_argument(
        "-f",
        "--filetype",
        type=str,
        choices=[
            'FASTQ',
            'SAM'],
        default='FASTQ',
        help="File type of the provided input data. Can be 'FASTQ' or 'SAM'."
        "Default is 'FASTQ'.")

    # Parse the arguments
    args = parser.parse_args()

    if args.filetype == 'FASTQ':
        while True:
            lines = []
            for _ in range(4):
                lines.append(args.input_file.readline())
            if (len(lines[0]) == 0):
                break
            lines[0] = lines[0].split(' ')[0]
            if lines[0][-1] != "\n":
                lines[0] += "\n"
            for _ in range(4):
                args.output_file.write(lines[_])

    elif args.filetype == 'SAM':
        for l in args.input_file:
            if l.startswith('@'):
                args.output_file.write(l)
            else:
                line = l.split('\t')
                line[0] = line[0].split(' ')[0]
                args.output_file.write("\t".join(line))

    args.input_file.close()
    args.output_file.close()


if __name__ == '__main__':
    main()
