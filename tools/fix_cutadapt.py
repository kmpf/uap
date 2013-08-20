#!/usr/bin/env python
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='This script reads one or two fastq files, line-by-line '+
        'and removes reads or read pairs if at least one read is empty.',
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("in_file1",
                        help="This file or fifo contains the " +
                        "first set of reads.")

    parser.add_argument("--in_file2",
                        dest="in_file2",
                        default=None,
                        help="This file or fifo contains the " +
                        "second set of reads.")

    parser.add_argument("out_file1",
                        help="To this file or fifo are the reads from " +
                        "in-file1 written without empty lines.")

    parser.add_argument("--out_file2",
                        dest="out_file2",
                        default=None,
                        help="This file or fifo contains the " +
                        "second set of reads.")

    args = parser.parse_args()        

    fin1 = open(args.in_file1, 'r')
    fin2 = args.in_file2
    fout1 = open(args.out_file1, 'w')
    fout2 = args.out_file2

    if args.in_file2 is not None:
        fin2 = open(args.in_file2, 'r')

    if args.out_file2 is not None:
        fout2 = open(args.out_file2, 'w')

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
            sys.stderr.write("Read %d entries, wrote %d entries.\n" % (rcount, wcount))

    sys.stderr.write("Read %d entries, wrote %d entries.\n" % (rcount, wcount))
        
if __name__ == '__main__':
    main()
