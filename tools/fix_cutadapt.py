#!/usr/bin/env python
import sys

def main():
    if len(sys.argv) != 5:
        print("Usage: fix_cutadapt.py [R1-in.fastq] [R2-in.fastq] [R1-out.fastq] [R2-out.fastq]")
        print("All files are uncompressed FASTQ files.")
        exit(1)

    fin1 = open(sys.argv[1], 'r')
    fin2 = open(sys.argv[2], 'r')
    fout1 = open(sys.argv[3], 'w')
    fout2 = open(sys.argv[4], 'w')

    rcount = 0
    wcount = 0
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

    sys.stderr.write("Read %d entries, wrote %d entries.\n" % (rcount, wcount))
        
if __name__ == '__main__':
    main()
