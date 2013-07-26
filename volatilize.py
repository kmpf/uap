#!./python_env/bin/python

import sys
sys.path.append('./include')
import argparse
import pipeline
import string
import yaml

'''
This script 
'''

parser = argparse.ArgumentParser(
    description="This script .",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("--details",
                    dest="details",
                    action="store_true",
                    help="Displays detailed information about ")

parser.add_argument("--srsly",
                    dest="srsly",
                    action="store_true",
                    help="")

parser.add_argument("--even-if-dirty",
                    dest="even_if_dirty",
                    action="store_true",
                    help="Must be set if the local git repository " +
                    "contains uncommited changes. Otherwise the pipeline " +
                    "will not start.")

args = parser.parse_args()

def main():
    p = pipeline.Pipeline()
    p.check_volatile_files(details = ('--details' in sys.argv), srsly = ('--srsly' in sys.argv))
        
if __name__ == '__main__':
    main()
