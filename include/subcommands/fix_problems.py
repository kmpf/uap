#!/usr/bin/env python

import sys
import os
import logging
import string
import yaml

import pipeline
from uaperrors import UAPError
'''
This script checks if anything went wrong with some tasks of the pipeline. It can
display more information if the '--details' is given and it can solve the problem
if '--srsly' is given.
'''

logger = logging.getLogger("uap_logger")

def main(args):
    p = pipeline.Pipeline(arguments=args)
 
    if args.first_error:
        report_script = p.get_cluster_command('last_error')
        cmd = [os.path.join(args.uap_path, report_script)]
        if report_script == '':
            raise UAPError('A first error report script is not implemented '
                    'for the cluster "%s".' % p.get_cluster_type())
        elif not os.path.exists(cmd[0]):
            raise UAPError('The configured report script "%s" does not exist.'
                    % cmd[0])
        ids = p.get_cluster_job_ids()
        if len(ids)==0:
            raise UAPError('No jobs found.')
        cmd.extend(ids)
        try:
            process = subprocess.call(cmd)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise UAPError('The configured first error report script "%s" '
                        'closed with an error.' % report_script)
            else:
                raise e
    else:
        p.check_ping_files(print_more_warnings = True,
                print_details = args.details, fix_problems = args.srsly)

#if __name__ == '__main__':
#    main()
