#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import yaml

p = pipeline.Pipeline()
p.print_tasks()
