#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline

p = pipeline.Pipeline()

p.print_tasks()
