#!/usr/bin/env python
import sys
import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
CHECK_QNAMES = False

class Filter:
    def __init__(self, source, destination):
        self.current_qname = None
        self.lines = list()
        self.flags_or = 0
        self.source = source
        self.destination = destination
        self.rcount = 0
        self.wcount = 0
        self.discard_40_80 = 0
        self.discard_gt_2 = 0
        if CHECK_QNAMES:
            self.seen_qnames = set()

    def handle(self):
        if len(self.lines) > 0:
            if CHECK_QNAMES and self.current_qname in self.seen_qnames:
                raise StandardError("Input is not sorted by qname. KA-BLAM.")

            if len(self.lines) <= 2:
                skip = False
                if len(self.lines) == 2:
                    if not ((self.flags_or & (0x40 | 0x80)) == (0x40 | 0x80)):
                        self.discard_40_80 += 1
                        skip = True
                if not skip:
                    if CHECK_QNAMES:
                        self.seen_qnames.add(self.current_qname)

                    for line in self.lines:
                        fixed_line = line
                        if len(self.lines) == 1:
                            line_array = line.split("\t")
                            flags = int(line_array[1])
                            if (flags & 0x8) == 0:
                                flags |= 0x8
                            line_array[1] = str(flags)
                            fixed_line = "\t".join(line_array)

                        self.destination.write(fixed_line)
                    self.wcount += len(self.lines)
            else:
                self.discard_gt_2 += 1

        self.current_qname = None
        self.lines = list()
        self.flags_or = 0

    def run(self):
        for line in self.source:
            if line[0] == '@':
                self.destination.write(line)
                continue

            self.rcount += 1

            line_array = line.split('\t')
            flags = int(line_array[1])

            if (flags & 0x100) != 0:
                continue

            qname = line_array[0]

            if qname != self.current_qname:
                self.handle()

            self.current_qname = qname
            self.lines.append(line)
            self.flags_or |= flags

        self.handle()
        self.wpercent = self.wcount * 100.0 / self.rcount
        keys = ['rcount', 'wcount', 'wpercent', 'discard_gt_2', 'discard_40_80']
        sys.stderr.write("#%s\n" % "\t".join(keys))
        sys.stderr.write("%s\n" % "\t".join([str(getattr(self, k)) for k in keys]))

def main():
    filter = Filter(sys.stdin, sys.stdout)
    filter.run()

if __name__ == '__main__':
    main()
