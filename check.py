#!./python_env/bin/python

# This script loads the configuration and prints some information about it.

import sys
sys.path.append('./include')
import pipeline
import yaml

p = pipeline.Pipeline()
print("Configuration is good, no errors were found.")

if len(sys.argv) > 1 and sys.argv[1] == '--detailed':
    print(yaml.dump(p.all_samples))
    exit(0)

print(p)
