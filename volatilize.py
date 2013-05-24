#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import string
import yaml

def main():
    p = pipeline.Pipeline()
    p.check_volatile_files(srsly = ('--srsly' in sys.argv))
        
if __name__ == '__main__':
    main()
