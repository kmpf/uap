#!./python_env/bin/python

# This script loads the configuration and prints some information about it.

import sys
sys.path.append('./include')
import pipeline

p = pipeline.Pipeline()
print >> sys.stderr, p
