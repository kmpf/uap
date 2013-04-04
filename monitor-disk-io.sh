#!/bin/bash
(strace -o /dev/stdout -f $*) | grep "write(" | grep -v "write(1" | grep -v "write(2"