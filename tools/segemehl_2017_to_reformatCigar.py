#!/bin/bash
"exec" "`dirname $0`/../python_env/bin/python" "$0" "$@"

# ^^^
# the cmd above ensures that the correct python environment is 
# selected to execute this script.
# The correct environment is the one belonging to uap, since all 
# neccessary python modules are installed there.


# filename: segemehl_2017_reformatCigar.py
# author: Jana Hertel
# date: 2017/05/31
# version: 1.0
# description: Format the .sam output produced by segemehl/export.13 module
#              into a .sam file that can be processed by htseq-count



import argparse
import sys
import re

parser = argparse.ArgumentParser(description='Python script to format the '
                                 '.sam output produced by segemehl/export.13 '
                                 'module into a .sam file that can be processed '
                                 'by htseq-count. The program writes to STDOUT.')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--in-file', dest='my_sam_in', required=True, 
                    type=argparse.FileType('r'),
                    help='Segemehl/export.13 .sam file that should be reformatted. '
                    'When piping in with samtools than just use - as argument.')

args = parser.parse_args()


################################################################################
# my_range(start, end, step)
#
# This function creates a range with an user defined step to walk through.
# returns: the respective new start values

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

################################################################################

for line in args.my_sam_in:

    columns = line.strip().split('\t')

    # don't process header lines
    # print to stdout as they are
    if(columns[0][:1]=="@"):
        for fields in columns:
            sys.stdout.write(fields)
            if(not fields==columns[len(columns)-1]):
                sys.stdout.write("\t")
        sys.stdout.write("\n")
    else: # all non-header lines
        # select cigar string
        cigar=columns[5]
        x = re.split('(\D)', cigar)

        # split cigar string and sum up values for '=' and 'X' (match and mismatch)
        # leave values as they are for 'I','D' and 'N' (del, insertion, split)
        M = 0
        cigar_new = ''
        for j in my_range(1, len(x)-1, 2):
            # match or mismatch
            if x[j]=='=' or x[j]=='X' or x[j]=='M':
                M = M + int(x[j-1])
            else: # del or ins
                cigar_new += str(M) + "M"
                M = 0
                cigar_new += x[j-1] + x[j]

        if M > 0:
            cigar_new += str(M) + "M"

        # print the sam line with the new cigar string to stdout
        for k in range(0, 4):
            sys.stdout.write("%s\t" % columns[k])

        sys.stdout.write("%s\t" % cigar_new)

        for k in range(6, len(columns)):
            sys.stdout.write("%s" % columns[k])
            if(not k==len(columns)):
                sys.stdout.write("\t")

        sys.stdout.write("\n")
