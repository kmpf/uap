#!/usr/bin/env python
import sys
import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
exec(compile(open(activate_this_file).read(), activate_this_file, 'exec'), dict(__file__=activate_this_file))
import random
import argparse
import gzip

def read_args():
    parser = argparse.ArgumentParser(description='randomly sample from [s]ingle or [p]aired end fastq file')
    parser.add_argument('readtype', type=str, metavar='type', nargs=1, choices=['single', 'paired']
                        , help='single or paired, default = None')
    parser.add_argument('--sample_size', '-n', type=int, nargs=1,              help='Number of reads to sample')
    parser.add_argument('--infiles', '-i', nargs='+',             help='Fastq input')
    parser.add_argument('--outfiles', '-o', nargs='+',             help='Fastq output')
    parser.add_argument('--read-gz',       nargs = '?', const=True, help='input is gzipped: infile.fastq.gz')
    parser.add_argument('--write-gz',       nargs='?', const=True, help='write gzipped output: outfile.fastq.gz')
    parser.add_argument('--fasta',          nargs='?', default='fastq', const='fasta', help='fasta or fastq file default is set to fastq')



    args = parser.parse_args()


    if args.fasta == 'fasta':
        args.divideby = 2
    else:
        args.divideby = 4

    #some checking
    if args.readtype[0] == 'single':
        if not (len(args.infiles)  == 1 and len(args.outfiles)  == 1):
            sys.stderr.write( "single: only 1 input and 1 output file allowed\n")
            exit(1)



    if args.readtype[0] == 'paired':
        if not (len(args.infiles)  == 2 and  len(args.outfiles)  == 2):
            sys.stderr.write( "paired: only 2 input and 2 output files allowed\n")
            exit(1)


    mylist = args.infiles + args.outfiles
    if (len(mylist)!=len(set(mylist))):
        sys.stderr.write( "duplicate filename\n")
        exit(1)

    return args


def provide_FH(args):
    my_range = {'single':1, 'paired':2}
    files = {'in':[], 'out':[]}

    for i in range(my_range[args.readtype[0]]):
        infilename = args.infiles[i]
        outfilename = args.outfiles[i]


        if args.read_gz:
            f = gzip.open(infilename, 'rb')
        else:
            f = open(infilename, 'r')
        files['in'].append(f)

        if args.write_gz:
            f = gzip.open(outfilename, 'wb')
        else:
            f = open(outfilename, 'w')
        files['out'].append(f)
    return files


def main(args):
    files = provide_FH(args)

    #read number of lines

    lines = 0
    buf_size = 1024 * 1024
    read_f = files['in'][0].read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    files['in'][0].seek(0, 0)



    records =  int(lines/args.divideby)

    if records < args.sample_size[0]:
        sys.stderr.write( "sample bigger than records aboarting\n")
        exit(1)

    rand_records = sorted(random.sample(range(records), args.sample_size[0]))



    rec_no = 0

    for rr in rand_records:
        while rec_no < rr:
             rec_no += 1
             for i in range(args.divideby): files['in'][0].readline()
             if args.readtype[0] == 'paired':
                 for i in range(args.divideby): files['in'][1].readline()

        for i in range(args.divideby):
            files['out'][0].write(files['in'][0].readline())
            if args.readtype[0] == 'paired':
                files['out'][1].write(files['in'][1].readline())

        rec_no += 1


if __name__ == '__main__':
    args = read_args()
    main(args)

