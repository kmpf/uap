#!/usr/bin/env python

import argparse
import os
import sys
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
    print("Reading .npz-file: %s" % f)
    with np.load(f, mmap_mode='r') as arr:
        labels.append(str(arr['labels'][0]))
        if not isinstance(ma, np.ndarray):
            ma = arr['matrix']
            print("Length arr['matrix']: %s" % len(arr['matrix']))
            print("Length ma: %s" % len(ma))
        else:
            readMatrix = arr['matrix']
            if len(readMatrix) < len(ma):
                print("%s" % (len(ma) - len(readMatrix)))
                a = np.array([[0]] * (len(ma) - len(readMatrix)))
                print("%s %s" % (a.shape, readMatrix.shape))
                readMatrix = np.vstack((readMatrix, a))
                print("%s %s" % (a.shape, readMatrix.shape))
            ma = np.concatenate((ma, readMatrix), axis=1)

        print("Matrix shape: %s" % str(ma.shape))
        print("Labels: %s" % ", ".join(labels))

    print(ma.size)

outFile = args.output
n = 0
while os.path.exists(outFile):
    n += 1
    head, tail = os.path.split(outFile)
    tail = "%s-%s" % (n, tail)
    outFile = os.path.join(head, tail)

print("Output file: %s" % outFile)
np.savez_compressed(outFile, matrix=ma, labels=labels)

    #print(arr)
    
    # check array dimensions if they are equal merge them by column    

    # np.load()
