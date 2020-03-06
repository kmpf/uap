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

    if p.task_wish_list:
        task_list = p.task_wish_list
    else:
        task_list = p.all_tasks_topologically_sorted

    # try to generate reports for all tasks
    for task in task_list:
        basic_task_state = task.get_task_state()
        if basic_task_state == p.states.FINISHED:
            try:
                task.generate_report()
            except:
                logger.info("Task %s did not produce")
        else:
            sys.stderr.write("Skipping %s because it's not finished yet.\n" % task)
