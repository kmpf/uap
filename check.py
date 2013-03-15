#!./python_env/bin/python

# This script loads the configuration and prints some information about it.

import sys
sys.path.append('./include')
import pipeline
import yaml

p = pipeline.Pipeline()

if len(sys.argv) > 1 and sys.argv[1] == '--detailed':
    print(yaml.dump(p.all_samples))
else:
    print(p)
    print("Configuration is good, no errors were found.")
    