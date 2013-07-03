#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import string
import yaml

def main():
    p = pipeline.Pipeline()
    p.check_ping_files(print_more_warnings = True, print_details = ('--details' in sys.argv), fix_problems = ('--srsly' in sys.argv))
        
if __name__ == '__main__':
    main()