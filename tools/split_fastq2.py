import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
exec(compile(open(activate_this_file).read(), activate_this_file, 'exec'), dict(__file__=activate_this_file))
import argparse
from Bio import SeqIO

'''

'''

def read_args():
    parser = argparse.ArgumentParser(
        description='generates splits of a fastq file'
    )
    parser.add_argument(
        '--infile',
        '-i',
        nargs=1,
        help='Fastq input'
    )
    parser.add_argument(
        '--read_count',
        '-n',
        type=int,
        nargs=1,
        help='Number of reads per output file'
    )
    parser.add_argument(
        '--outpath',
        '-o',
        nargs='?',
        help='output path'
    )
    parser.add_argument(
        '--out_file_pattern',
        '-p',
        nargs='?',
        help='file name pattern for output files, a counter will be added automatically at the end (<pattern>_1.fastq)'
    )
    parser.add_argument(
        '--sindex',
        '-s',
        nargs='?',
        help='index of sample to write to file'
    )
    parser.add_argument(
        '--mate',
        '-m',
        nargs='?',
        help='mate1: r1, mate2: r2'
    )

    args = parser.parse_args()
    # TODO: catch missing args

    return args

def main(args):
    # TODO: Exception if its not a fastq
    infile = args.infile[0]
    read_count = args.read_count[0]
    out_file_pattern = args.out_file_pattern if args.out_file_pattern is not None else 'split'

    # TODO: check for valid and existing path
    outpath = args.outpath if args.outpath is not None else ''
    if outpath != '' and outpath[-1] != '/':
        outpath += '/'

    sample_index = args.sindex
    mate = args.mate
    chunksize = read_count * 4

    start = (int(sample_index) - 1) * int(chunksize)
    with open(infile, 'r') as fin:
        tmp_filename = outpath + "%s_%i_%s.fastq" % (out_file_pattern, int(sample_index), mate)
        fout = open(tmp_filename, "wb")
        for i, line in enumerate(fin):
            if i >= start:
                fout.write(line)
                if (i + 1) % chunksize == 0:
                    fout.close()
                    break
        fout.close()

    # if last filename is empty delete
    filesize = os.path.getsize(tmp_filename)
    if filesize == 0:
        os.remove(tmp_filename)

if __name__ == '__main__':
    args = read_args()
    main(args)