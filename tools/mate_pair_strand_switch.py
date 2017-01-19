#!/usr/bin/env python
import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import sys
import argparse
import fileinput
import re

def main():
    parser = argparse.ArgumentParser(
        description='This script reads a bed file created by \'bedtools bamtobed\' ' +
        'line-by-line and reverts the strand information for every mate pair read.',
        prog = 'mate_pair_strand_switch.py',
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 0.01 alpha')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin, help ="Infile default reads from stdin")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stdout, help ="Outfile default writes to stdout")


    args = parser.parse_args()

    for bed_line in args.infile:
        bed_line = bed_line.rstrip()
        columns = bed_line.split('\t')
        if re.search('2$', columns[3]):
            columns[5] = '-' if columns[5] == '+' else '+'
        sys.stdout.write('\t'.join(columns) + '\n')    

if __name__ == '__main__':
    main()
