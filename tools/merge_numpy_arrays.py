#!/usr/bin/env python

import argparse

import numpy as np
from scipy import sparse

parser = argparse.ArgumentParser(
    add_help=True,
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    "files",
    nargs="*",
    default=list(),
    type=str,
    help="All .npz files which are merged"
)

parser.add_argument(
    "output",
    type=str,
    help="Output file basename"
)

args = parser.parse_args()

ma = None
labels = list()
for f in args.files:
    with np.load(f, mmap_mode='r') as arr:
        labels.append(str(arr['labels'][0]))
        if not isinstance(ma, np.ndarray):
            ma = arr['matrix']
        else:
            ma = np.concatenate((ma, arr['matrix']), axis=1)

        print("Matrix shape: %s" % str(ma.shape))
        print("Labels: %s" % ", ".join(labels))

    print(ma.size)

np.savez_compressed(args.output, matrix=ma, labels=labels)

    #print(arr)
    
    # check array dimensions if they are equal merge them by column    

    # np.load()
