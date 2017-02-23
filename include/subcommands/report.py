#!/usr/bin/env python

import sys
import copy
import logging
import os
import signal
import socket
import yaml

import misc
import pipeline
import process_pool

logger = logging.getLogger("uap_logger")

def main(args):
    p = pipeline.Pipeline(arguments=args)
    def handle_signal(signum, frame):
        print("Catching %s!" % process_pool.ProcessPool.SIGNAL_NAMES[signum])
        p.caught_signal = signum
        process_pool.ProcessPool.kill()

    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)

    task_list = p.all_tasks_topologically_sorted

    if len(args.run) >= 1:
        # execute the specified tasks
        task_list = list()
        for task_id in args.run:
            if '/' in task_id:
                task = p.task_for_task_id[task_id]
                task_list.append(task)
            else:
                for task in p.all_tasks_topologically_sorted:
                    if str(task)[0:len(task_id)] == task_id:
                        task_list.append(task)

    # try to generate reports for all tasks
    for task in task_list:
        basic_task_state = task.get_task_state_basic()
        if basic_task_state == p.states.FINISHED:
            try:
                task.generate_report()
            except:
                logger.info("Task %s did not produce")
        else:
            sys.stderr.write("Skipping %s because it's not finished yet.\n" % task)
