#!/usr/bin/env python

import sys
import os
import logging
import string
import yaml
import subprocess
from tqdm import tqdm
from datetime import datetime

import pipeline
from uaperrors import UAPError
'''
This script checks if anything went wrong with some tasks of the pipeline. It can
display more information if the '--details' is given and it can solve the problem
if '--srsly' is given.
'''

logger = logging.getLogger("uap_logger")

def main(args):
    args.no_tool_checks = True
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
    elif args.file_modification_date:
        if not args.srsly:
            print('This is a dry run and no files will be changed. '
                  'Pass --srsly to apply changes.')
        tasks = p.all_tasks_topologically_sorted
        question = 'This will change the modification dates of all output ' \
                   'files with valid sha256sum. Do you want to proceed?'
        split_at = ' sha256sum is correct and modification date after '
        split_at_v = ' is volatilized and modification date after '
        changes = list()
        task_iter = tqdm(tasks, desc='tasks')
        try:
            for task in task_iter:
                for bad_file in task.get_run().file_changes(do_hash=True):
                    bad_file = bad_file.split(split_at)
                    if len(bad_file) == 1:
                        bad_file = bad_file[0].split(split_at_v)
                    if len(bad_file) == 1:
                        continue
                    date = task.get_run().written_anno_data()['end_time']
                    changes.append('%s set for %s' % (date, bad_file[0]))
                    if args.srsly:
                        os.utime(bad_file[0], (now(), timestamp(date)))
        except:
            task_iter.close()
            raise
        if changes:
            print('#### Changes ####')
        else:
            print('Nothing chages.')
        for line in changes:
            print(line)

    else:
        p.check_ping_files(print_more_warnings = True,
                print_details = args.details, fix_problems = args.srsly)

def timestamp(date):
    return (date - datetime.fromtimestamp(0)).total_seconds()

def now():
    return timestamp(datetime.now())

#if __name__ == '__main__':
#    main()
