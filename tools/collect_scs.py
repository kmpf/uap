#! /usr/bin/env python

# collect_metrics_seqXTo2d.py
# Author: Sven Holger Pupple
# This script collects metrics from successive bowtie2
# and count_rRNA steps.

import sys
import re
import datetime
import argparse
import yaml


def read_arguments():
    parser = argparse.ArgumentParser(description="scs succesive screen")
    parser.add_argument('--infiles', nargs='+', type=argparse.FileType('r'),
                        default=sys.stdin, help="first alignment summary")
    parser.add_argument('--outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help="Outfile: default=stdout")
    parser.add_argument('--types', nargs='+', type=str,
                        default=None, help="Outfile: default=stdout")
    parser.add_argument('--rrna-aln-pos', nargs='?', type=str, default=None,
                        help="rRNA counts")
    choices = ['unstranded', 'firststranded', 'secondstranded']
    parser.add_argument('--library-type', choices=choices,
                        default=None, help="See tophat manual.")
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('--name', default=None, type=str, help="set id")
    return parser.parse_args()


def process_bowtie2(f):
    res = dict()
    keys = dict()

    p_zero = re.compile(r'0\stimes')
    p_once = re.compile(r'exactly 1 time')
    p_multiple = re.compile(r'>1 times')

# single
# 9202262 reads; of these:
#   9202262 (100.00%) were unpaired; of these:
#     859754 (9.34%) aligned 0 times
#     1233838 (13.41%) aligned exactly 1 time
#     7108670 (77.25%) aligned >1 times
# 90.66% overall alignment rate
# paired
# 3826398 reads; of these:
#   3826398 (100.00%) were paired; of these:
#     3826396 (100.00%) aligned concordantly 0 times
#     0 (0.00%) aligned concordantly exactly 1 time
#     2 (0.00%) aligned concordantly >1 times
#     ----

    for line in f:
        if not line.strip():
            continue
        line = line.rstrip('\n')
        line = line.lstrip(' ')
        tmp = line.split(' ')

        count = tmp[0]
        other = tmp[1]

        if count == '----':
            break
        if other == 'reads;':
            pass

        if tmp[2] == 'were':
            res['input'] = int(count)
            keys['input'] = 1

        if p_zero.search(line):
            keys['unaligned'] = 1
            keys['p_unaligned'] = 1
            keys['overall'] = 1
            res['overall'] = res['input'] - int(count)
            res['unaligned'] = count
            res['p_unaligned'] = other.split('%')[0][1:]

        if p_once.search(line):
            keys['unique'] = 1
            keys['p_unique'] = 1
            res['unique'] = count
            res['p_unique'] = other.split('%')[0][1:]

        if p_multiple.search(line):
            keys['multiple'] = 1
            keys['p_multiple'] = 1
            res['multiple'] = count
            res['p_multiple'] = other.split('%')[0][1:]

        if other == 'overall':
            keys['p_overall'] = 1
            res['p_overall'] = count.split('%')[0]

    ress = dict()
    ress['unique'] = res['unique']
    ress['multiple'] = res['multiple']
    ress['unaligned'] = res['unaligned']

    return ress


def process_count_rRNA(args, f):

    # Illumina Trueseq firststranded
    # by definition
    # first read reverse -> comes from +     samflag 16
    # first read not reverse -> comes from - samflag 0
    # second read reverse -> comes from -    samflag 16
    # second not read reverse -> comes from + samflag 0

    # truseq the reverse

    res = dict()
    res2 = dict()
    res['sense'] = 0
    res['antisense'] = 0
    res['unstranded'] = 0

    for line in f:
        if not line.strip():
            continue
        line = line.rstrip('\n')
        line = line.lstrip(' ')
        (count, flag, ref) = line.split()

        if flag == '0' or flag == '16':
            if args.library_type == 'unstranded':
                res['unstranded'] += int(count)

            if args.library_type == 'firststranded':
                if flag == '0':
                    t = "-".join(['antisense', ref])
                    res2[t] = int(count)
                    res['antisense'] += int(count)
                else:
                    t = "-".join(['sense', ref])
                    res2[t] = int(count)
                    res['sense'] += int(count)

            if args.library_type == 'secondstranded':
                if flag == '16':
                    t = "-".join(['antisense', ref])
                    res2[t] = int(count)
                    res['antisense'] += int(count)
                else:
                    t = "-".join(['sense', ref])
                    res2[t] = int(count)
                    res['sense'] += int(count)
        else:
            print(count, flag, ref)
            raise Exception("Other flag than 0 or 16 encountered.")

    return res, res2


def main(args):
    outd = dict()
    outd['id'] = args.name
    outd['programm'] = "scs_collect"
    outd['date'] = datetime.datetime.now()

    for n, f in enumerate(args.infiles):
        outd[n] = dict()
        outd[n]['type'] = args.types[n]

        if str(n + 1) == args.rrna_aln_pos:
            (outd[n]['values'], outd['rnacomposition']) = \
                    process_count_rRNA(args, f)
        else:
            outd[n]['values'] = process_bowtie2(f)
            # remove unaligned for all except last entry
            if n + 1 != len(args.infiles):
                del outd[n]['values']['unaligned']

    args.outfile.write(yaml.dump(outd, default_flow_style=False,
                       default_style=''))


__version__ = '0.001'

if __name__ == '__main__':
    args = read_arguments()
    main(args)
