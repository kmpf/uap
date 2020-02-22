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
        if args.srsly and not ask_user(question):
            return
        split_at = ' sha256sum is correct and was changed after '
        changes = list()
        for task in tqdm(tasks, desc='tasks'):
            for bad_file in task.get_run().file_changes(do_hash=True):
                if not isinstance(bad_file, str):
                    break
                bad_file = bad_file.split(split_at)
                if len(bad_file) == 1:
                    continue
                date = task.get_run().written_anno_data()['end_time']
                changes.append('%s set for %s' % (date, bad_file[0]))
                if args.srsly:
                    os.utime(bad_file[0], (now(), timestamp(date)))
        if changes:
            print('#### Changes ####')
        else:
            print('Nothing chages.')
        for line in changes:
            print(line)

    else:
        p.check_ping_files(print_more_warnings = True,
                print_details = args.details, fix_problems = args.srsly)

def ask_user(question):
    check = str(raw_input("%s (y/n): " % question)).lower().strip()
    try:
        if check[0] == 'y':
            return True
        elif check[0] == 'n':
            return False
        else:
            print('Invalid Input')
            return ask_user(question)
    except Exception as error:
        print("Please enter valid inputs")
        print(error)
        return ask_user(question)

def timestamp(date):
    return (date - datetime.fromtimestamp(0)).total_seconds()

def now():
    return timestamp(datetime.now())

#if __name__ == '__main__':
#    main()
