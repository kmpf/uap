#!/usr/bin/env python
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

    args = parser.parse_args()

    for bed_line in fileinput.input():
        bed_line = bed_line.rstrip()
        columns = bed_line.split('\t')
        if re.search('2$', columns[3]):
            columns[5] = '-' if columns[5] == '+' else '+'
        print('\t'.join(columns))    

if __name__ == '__main__':
    main()
