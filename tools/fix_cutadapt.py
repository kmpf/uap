#!/bin/bash
import argparse
import sys
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
    parser = argparse.ArgumentParser(
        description='This script reads one or two fastq files, line-by-line ' +
        'and removes reads or read pairs if at least one read is empty.',
        prog='fix_cutadapt',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 0.01')

    parser.add_argument("r1_in",
                        help="This file or fifo contains the " +
                        "first set of reads.")

    parser.add_argument("--R2-in",
                        dest="r2_in",
                        default=None,
                        help="This file or fifo contains the " +
                        "second set of reads.")

    parser.add_argument("r1_out",
                        help="To this file or fifo are the reads from " +
                        "in-file1 written without empty lines.")

    parser.add_argument("--R2-out",
                        dest="r2_out",
                        default=None,
                        help="This file or fifo contains the " +
                        "second set of reads.")

    args = parser.parse_args()

    fin1 = open(args.r1_in, 'r')
    fin2 = None
    fout1 = open(args.r1_out, 'w')
    fout2 = None

    if args.r2_in is not None:
        fin2 = open(args.r2_in, 'r')

    if args.r2_out is not None:
        fout2 = open(args.r2_out, 'w')

    rcount = 0
    wcount = 0
    if fin2 is not None and fout2 is not None:
        while True:
            lines1 = []
            lines2 = []
            for _ in range(4):
                lines1.append(fin1.readline())
                lines2.append(fin2.readline())
            if (len(lines1[0]) == 0):
                break
            if not (lines1[1] == "\n" or lines2[1] == "\n"):
                for _ in range(4):
                    fout1.write(lines1[_])
                    fout2.write(lines2[_])
                wcount += 1
            rcount += 1

        fin1.close()
        fin2.close()
        fout1.close()
        fout2.close()

    else:
        while True:
            lines1 = []
            for _ in range(4):
                lines1.append(fin1.readline())
            if (len(lines1[0]) == 0):
                break
            if not (lines1[1] == "\n"):
                for _ in range(4):
                    fout1.write(lines1[_])
                wcount += 1
            rcount += 1
            sys.stderr.write(
                "Read %d entries, wrote %d entries.\n" %
                (rcount, wcount))

    sys.stderr.write("Read %d entries, wrote %d entries.\n" % (rcount, wcount))


if __name__ == '__main__':
    main()
