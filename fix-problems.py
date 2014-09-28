#!./python_env/bin/python

import sys
sys.path.append('./include')
import argparse
import pipeline
import string
import yaml

'''
This script checks if anything went wrong with some tasks of the pipeline. It can  
display more information if the '--details' is given and it can solve the problem
if '--srsly' is given.
'''

parser = argparse.ArgumentParser(
    description="This script displays by default information about all tasks " +
                "of the pipeline as configured in 'config.yaml'. But the " +
                "displayed information can be narrowed down via command " +
                "line options.",
    formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument("--cluster",
                    dest="cluster",
                    type=str,
                    default="auto",
                    help="Specify the cluster type (sge, slurm), defaults to auto.")


parser.add_argument("--details",
                    dest="details",
                    action="store_true",
                    default=False,
                    help="Displays detailed information about problems that " +
                    "occured during the execution of the pipeline.")

parser.add_argument("--srsly",
                    dest="srsly",
                    action="store_true",
                    default=False,
                    help="Fixes the problems that the pipeline has encountered.")

parser.add_argument("--even-if-dirty",
                    dest="even_if_dirty",
                    action="store_true",
                    default=False,
                    help="Must be set if the local git repository " +
                    "contains uncommited changes. Otherwise the pipeline " +
                    "will not start.")

args = parser.parse_args()

def main():
    p = pipeline.Pipeline(arguments=args)
    p.check_ping_files(print_more_warnings = True, print_details = args.details, fix_problems = args.srsly)
        
if __name__ == '__main__':
    main()
