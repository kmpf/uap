#!/usr/bin/env python

import sys
import logging
import string
import yaml

from .. import pipeline

'''
This script 
'''

logger = logging.getLogger("uap_logger")

def main():
    p = pipeline.Pipeline()
    p.check_volatile_files(details = args.details, srsly = args.srsly)
        
if __name__ == '__main__':
    main()
