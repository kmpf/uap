#!/usr/bin/env python

import sys
import logging
import string
import yaml

import pipeline

'''
This script
'''

logger = logging.getLogger("uap_logger")


def main(args):
    p = pipeline.Pipeline(arguments=args)
    p.check_volatile_files(details=args.details, srsly=args.srsly)

