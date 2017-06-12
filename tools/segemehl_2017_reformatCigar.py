#!/bin/bash
"exec" "`dirname $0`/../python_env/bin/python" "$0" "$@"

# ^^^
# the cmd above ensures that the correct python environment is 
# selected to execute this script.
# The correct environment is the one belonging to uap, since all 
# neccessary python modules are installed there.


# filename: segemehl_2017_reformatCigar.py
# author: Jana Hertel
# date: 2017/06/07
# version: 1.0
# description: Reformat the cigar string such that htseq-count is able to process
#              the according SAM files. Consecutive values for 'ix', 'j=' and 'kM'
#              are summed up and replaced by nM with n being the sum of i, j and k.


import argparse
import sys
import re

parser = argparse.ArgumentParser(description='Python script to process a large file '
                                 'using multi-processing.')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--in-file', dest='my_file_in', required=True, 
                    type=argparse.FileType('r'),
                    help='A large file whose lines are independent from each other and '
                    'can be processed separately.')
parser.add_argument('--threads', dest='my_cores', default=1,
                    type=int,
                    help='Number of CPUs 2B used.')

args = parser.parse_args()

################################################################################
# my_range(start, end, step)
#
# This function creates a range with a user defined step to walk through.
# returns: the respective new start values

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

################################################################################

################################################################################
# process_line(line)
#
# function that does something with the line:
# in this case:
# - split the line into columns by tab
# - print the columns separated by tab to stdout

def process_line(line):
#    sys.stdout.write("\n%s" % line)
    columns = line.strip().split('\t')
    cigar=columns[5]
    x = re.split('(\D)', cigar)

    # split cigar string and sum up consecutive values 
    # for '=' and 'X' (match and mismatch)
    # leave values as they are for 'I','D' and 'N' (del, insertion, split)
    M = 0
    cigar_new = ''
    for j in range(1, len(x)-1, 2):
        # match or mismatch
        if x[j]=='=' or x[j]=='X' or x[j]=='M':
            M = M + int(x[j-1])
        else: # del or ins
            if M > 0:
                cigar_new += str(M) + "M" # print the previous match/mismatch
                M = 0
            cigar_new += x[j-1] + x[j] # anything else but '=', 'X' or 'M'
            
    if M > 0:
        cigar_new += str(M) + "M"
        
    if cigar_new == "0M*":
        cigar_new = "*"

    # print the sam line with the new cigar string to stdout
    for k in range(0, 5):
        sys.stdout.write("%s\t" % columns[k])

    sys.stdout.write("%s\t" % cigar_new)

    for k in range(6, len(columns)):
        sys.stdout.write("%s" % columns[k])
        if(not k==len(columns)):
            sys.stdout.write("\t")

    sys.stdout.write("\n")

# END: process_line(line)
################################################################################





if __name__ == '__main__':
    
    for line in args.my_file_in:
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
            process_line(line)

    
