#!/usr/bin/env python
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
