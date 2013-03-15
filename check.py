#!./python_env/bin/python

# This script loads the configuration and prints some information about it.

import sys
sys.path.append('./include')
import pipeline
import yaml

p = pipeline.Pipeline()
print("")
print(p)
print("Detailed information:")
print(yaml.dump(p.all_samples))