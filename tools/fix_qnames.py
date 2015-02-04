#!/usr/bin/env python
import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import sys

def main():
    while True:
        lines = []
        for _ in range(4):
            lines.append(sys.stdin.readline())
        if (len(lines[0]) == 0):                                                                                                         
            break
        lines[0] = lines[0].split(' ')[0]
        if lines[0][-1] != "\n":
            lines[0] += "\n"
        for _ in range(4):
            sys.stdout.write(lines[_])

if __name__ == '__main__':
    main()
