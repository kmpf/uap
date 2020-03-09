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

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    From http://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

def main(args):
    # TODO: Exception if its not a fastq
    infile = args.infile[0]
    read_count = args.read_count[0]
    out_file_pattern = args.out_file_pattern if args.out_file_pattern is not None else 'group'

    # TODO: check for valid and existing path
    outpath = args.outpath if args.outpath is not None else ''
    if outpath != '' and outpath[-1] != '/':
        outpath += '/'

    sample_index = args.sindex
    mate = args.mate

    record_iter = SeqIO.parse(open(infile), "fastq")
    for i, batch in enumerate(batch_iterator(record_iter, read_count)):
        if int(sample_index) == (i+1):
            filename = outpath + "%s_%i_%s.fastq" % (out_file_pattern, int(sample_index), mate)
            with open(filename, "w") as handle:
                SeqIO.write(batch, handle, "fastq")

if __name__ == '__main__':
    args = read_args()
    main(args)