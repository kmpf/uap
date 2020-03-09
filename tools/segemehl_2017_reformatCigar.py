#!/bin/bash
"exec" "`dirname $0`/../python_env/bin/python" "$0" "$@"
#"exec" "python" "$0" "$@"

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
from multiprocessing import Pool
import itertools

parser = argparse.ArgumentParser(description='Python script to process a large file '
                                 'using multi-processing.')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--in-file', dest='my_file_in', required=True, 
                    type=argparse.FileType('r'),
                    help='A large file whose lines are independent from each other and '
                    'can be processed separately.')
parser.add_argument('--threads', dest='my_cores', default=1,
                    type=int,
                    help='Number of CPUs 2B used. Default: 1')
parser.add_argument('--blocksize', dest='my_bufsize', default=2,
                    type=int,
                    help='Size of buffer to read the input file (in MB). Default: 2')

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
# - returns the columns separated by tab 

def process_line(lines):

    newlines = list()

    c = 0
    for line in lines:
        c += 1

        columns = line.strip().split('\t')

        # don't process header lines
        if(columns[0][:1]=="@"):
            newlines.append(line.strip())
            continue
    
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
        new_line = ""
        for k in range(0, 5):
            new_line += "%s\t" % columns[k]

        new_line += "%s\t" % cigar_new

        for k in range(6, len(columns)):
            new_line += "%s" % columns[k]
            if(not k==len(columns)):
                new_line += "\t"

        newlines.append(new_line)
        
    return newlines

# END: process_line(line)
################################################################################




if __name__ == '__main__':

    # create my_cores -1 pools, 1 control + the remaining for processing the lines
    p = Pool(args.my_cores)

    a = list()

    eof_reached = False

    # bufsize needs to be provided in bytes
    # argument provided megabytes
    bufsize = args.my_bufsize * 1000000

    while not eof_reached:
        for i in range(args.my_cores - 1):
            linelist = args.my_file_in.readlines(bufsize)
            if len(linelist) == 0:
                eof_reached = True
            else:
                a.append(linelist) # ~ 2MB chunks
            
        l = p.map(process_line, a)

        for j in l:
            print('\n'.join(j))

        a[:] = [] # delete processed lines from the list


# this works in principle.. too much i/o   
#    for line in p.imap(process_line, args.my_file_in):
#        print line, # the coma prevents printing an additional new line


# idea for mp:
# read file in chunks of the size 1/args.my_cores
# --> each chunk in one process
