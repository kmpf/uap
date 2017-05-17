#!/usr/bin/env python

import sys
import logging
import glob
import os
import string
import yaml

import pipeline

'''
This script
'''

logger = logging.getLogger("uap_logger")

def main(args):
    p = pipeline.Pipeline(arguments=args)

    # Compile a list of all tasks
    task_list = list()
    # Only use tasks listed in args.run
    arguments = vars(args)

    if arguments['run_id'] and len(arguments['run_id']) >= 1:
        for task_id in arguments['run_id']:
            if '/' in task_id:
                task = p.task_for_task_id[task_id]
                task_list.append(task)
            else:
                for task in p.all_tasks_topologically_sorted:
                    if str(task)[0:len(task_id)] == task_id:
                        task_list.append(task)
    # or take all available tasks
    else:
        task_list = p.all_tasks_topologically_sorted

    for task in task_list:
        print(task)
        outdir = task.get_run().get_output_directory()
        anno_files = glob.glob(os.path.join(
            outdir, ".%s*.annotation.yaml" % task.get_run().get_run_id()
        ))

        yaml_files = {os.path.realpath(f) for f in anno_files \
                      if os.path.islink(f)}
        for y in yaml_files:
            print(y)


if __name__ == '__main__':
    main(args)
