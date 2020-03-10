#!/usr/bin/env python

import sys
import copy
import logging
import os
import signal
import socket
import yaml
from datetime import datetime

import misc
import pipeline
import process_pool
from uaperrors import UAPError

logger = logging.getLogger("uap_logger")


def main(args):
    p = pipeline.Pipeline(arguments=args)

    task = None
    def handle_signal(signum, frame):
        logger.warning("Catching %s!" %
                       process_pool.ProcessPool.SIGNAL_NAMES[signum])
        p.caught_signal = signum
        process_pool.ProcessPool.kill()
        if task:
            signame = process_pool.ProcessPool.SIGNAL_NAMES[signum]
            error = 'UAP stopped because it caught signal %d - %s' % \
                        (signum, signame)
            log_task_error(task, error, True)
    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)

    # execute all tasks
    finished_states = [p.states.FINISHED]
    if args.ignore:
        finished_states += [p.states.CHANGED]

    accepted_states = [p.states.BAD, p.states.READY, p.states.QUEUED,
                       p.states.VOLATILIZED]
    for task in p.get_task_with_list():
        task_state = task.get_task_state()
        if task_state in finished_states:
            task.move_ping_file()
            sys.stderr.write("Skipping %s because it's already %s.\n" %
                             (task, task_state))
        elif task_state == p.states.VOLATILIZED and not args.run:
            task.move_ping_file()
            sys.stderr.write("Skipping %s because it's already %s and not "
                             "specified as argument.\n" %
                             (task, task_state))
        elif task_state == p.states.CHANGED:
            if not args.force:
                task.move_ping_file()
                raise UAPError(
                    "Task %s has changed. "
                    "Run 'uap %s status --details' to see what changed or "
                    "'uap %s run-locally --force' to force overwrite "
                    "of the results." %
                    (task, args.config.name, args.config.name))
            else:
                check_parents_and_run(task, finished_states, args.debugging)
        elif task_state in accepted_states:
            check_parents_and_run(task, finished_states, args.debugging)
        else:
            task.move_ping_file()
            raise UAPError(
                "Unexpected task state for %s: %s\n"
                "Expected state to be 'READY'. Probably an upstream "
                "run crashed." %
                (task, task_state))


def check_parents_and_run(task, states, turn_bad):
    parents = task.get_parent_tasks()
    for parent_task in parents:
        parent_state = parent_task.get_task_state()
        if parent_state not in states:
            should = ' or '.join(states)
            error = "Cannot run %s because a parent job " \
                "%s is %s when it should be %s." % \
                (task, parent_task, parent_state, should)
            log_task_error(task, error, turn_bad)
    error = 'Task failed write a detailed annotation'
    try:
        task.run()
    except BaseException:
        error += ' and crashed with:\n%s' % \
                ''.join(traceback.format_exception(
                *sys.exc_info())[-2:]).strip()
    log_task_error(task, error, True)


def log_task_error(task, error, turn_bad):
    if turn_bad:
        run = task.get_run()
        run.get_step().start_time = datetime.now()
        run.get_step().end_time = datetime.now()
        if not run.annotation_written:
            run.write_annotation_file(error=error)
        task.move_ping_file()
    else:
        task.move_ping_file(bad_copy=False)
    raise UAPError(error)

