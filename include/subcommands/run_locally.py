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
from uaperrors import UAPError

logger = logging.getLogger("uap_logger")

def main(args):
    p = pipeline.Pipeline(arguments=args)

    original_term_handler = signal.getsignal(signal.SIGTERM)
    original_int_handler = signal.getsignal(signal.SIGINT)
    def handle_signal(signum, frame):
        logger.debug("Catching %s!" % process_pool.ProcessPool.SIGNAL_NAMES[signum])
        p.caught_signal = signum
        process_pool.ProcessPool.kill()
        if signum == signal.SIGTERM:
            original_term_handler(signum, frame)
        if signum == signal.SIGINT:
            original_int_handler(signum, frame)
        signal.SIG_DFL(signum, frame)
    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)

    #task_list = copy.deepcopy(p.all_tasks_topologically_sorted)
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
            
    # execute all tasks
    for task in task_list:
        task_state = task.get_task_state()
        if task_state == p.states.FINISHED:
            task.move_ping_file()
            sys.stderr.write("Skipping %s because it's already finished.\n" %
                             task)
        elif task_state == p.states.CHANGED:
            if args.ignore:
                task.move_ping_file()
                sys.stderr.write("Skipping %s because it's changes are "
                                 "ignored.\n" % task)
            elif not args.force:
                task.move_ping_file()
                raise UAPError("Task %s is finished but its config changed. "
                        "Run 'uap %s status --details' to see what changed or "
                        "'uap %s run-locally --force' to force overwrite "
                        "of the results." %
                        (task, args.config.name, args.config.name))
            else:
                task.run()
        elif task_state == p.states.BAD:
            if not args.force:
                task.move_ping_file()
                raise UAPError("Task %s is BAD. Resolve this problem with "
                        "'uap %s fix-problems' or fore an overwrite with "
                        "'uap %s run-locally --force'." %
                        (task, args.config.name, args.config.name))
            else:
                task.run()
        elif task_state == p.states.READY:
            task.run()
        elif task_state == p.states.QUEUED:
            task.run()
        else:
            task.move_ping_file()
            raise UAPError("Unexpected basic task state for %s: %s\n"
                         "Expected state to be 'READY'. Probably an upstream "
                         "run crashed." %
                         (task, task_state))

if __name__ == '__main__':
    try:
        main()
    finally:
        # make sure all child processes get terminated
        process_pool.ProcessPool.kill_all_child_processes()
