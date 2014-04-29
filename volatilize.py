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
                    default=False,
                    help="Displays detailed information about ")

parser.add_argument("--srsly",
                    dest="srsly",
                    action="store_true",
                    default=False,
                    help="")

parser.add_argument("--even-if-dirty",
                    dest="even_if_dirty",
                    action="store_true",
                    default=False,
                    help="Must be set if the local git repository " +
                    "contains uncommited changes. Otherwise the pipeline " +
                    "will not start.")

args = parser.parse_args()

def main():
    p = pipeline.Pipeline()
    p.check_volatile_files(details = args.details, srsly = args.srsly)
        
if __name__ == '__main__':
    main()
