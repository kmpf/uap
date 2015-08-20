#!/usr/bin/env python

import sys
import logging
import string
import yaml

from .. import pipeline
'''
This script checks if anything went wrong with some tasks of the pipeline. It can  
display more information if the '--details' is given and it can solve the problem
if '--srsly' is given.
'''

logger = logging.getLogger("uap_logger")

def main(args):
    p = pipeline.Pipeline(arguments=args)
    p.check_ping_files(print_more_warnings = True, print_details = args.details, fix_problems = args.srsly)
        
#if __name__ == '__main__':
#    main()
