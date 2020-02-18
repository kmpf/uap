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

    def handle_signal(signum, frame):
        logger.debug("Catching %s!" % process_pool.ProcessPool.SIGNAL_NAMES[signum])
        p.caught_signal = signum
        process_pool.ProcessPool.kill()
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
    if args.ignore:
        finished_states = [p.states.FINISHED, p.states.CHANGED]
    else:
        finished_states = [p.states.FINISHED]
    for task in task_list:
        task_state = task.get_task_state()
        if task_state in finished_states:
            task.move_ping_file()
            sys.stderr.write("Skipping %s because it's already %s.\n" %
                             (task, task_state))
        elif task_state == p.states.CHANGED:
            if not args.force:
                task.move_ping_file()
                raise UAPError("Task %s is finished but its config changed. "
                        "Run 'uap %s status --details' to see what changed or "
                        "'uap %s run-locally --force' to force overwrite "
                        "of the results." %
                        (task, args.config.name, args.config.name))
            else:
                check_parents_and_run(task, finished_states)
        elif task_state == p.states.BAD:
            if not args.force:
                raise UAPError("Task %s is BAD. Resolve this problem with "
                        "'uap %s fix-problems' or fore an overwrite with "
                        "'uap %s run-locally --force'." %
                        (task, args.config.name, args.config.name))
            else:
                check_parents_and_run(task, finished_states)
        elif task_state == p.states.READY:
            check_parents_and_run(task, finished_states)
        elif task_state == p.states.QUEUED:
            check_parents_and_run(task, finished_states)
        else:
            task.move_ping_file()
            raise UAPError("Unexpected basic task state for %s: %s\n"
                         "Expected state to be 'READY'. Probably an upstream "
                         "run crashed." %
                         (task, task_state))

def check_parents_and_run(task, states):
    parents = task.get_parent_tasks()
    for parent_task in parents:
        parent_state = parent_task.get_task_state()
        if parent_state not in states:
            raise UAPError("Cannot run %s because a parent job "
                           "%s is %s when it should be in %s." %
                           (task, parent_task, parent_state, states))
    task.run()

if __name__ == '__main__':
    try:
        main()
    finally:
        # make sure all child processes get terminated
        process_pool.ProcessPool.kill_all_child_processes()
