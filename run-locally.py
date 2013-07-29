#!./python_env/bin/python

import sys
sys.path.append('./include')
import argparse
import copy
import misc
import os
import pipeline
import process_pool
import signal
import socket
import yaml

parser = argparse.ArgumentParser(
    description="This script starts the 'rnaseq-pipeline' on the local machine. " +
                "It can be used to start:\n" +
                " * all tasks of the pipeline as configured in 'config.yaml'\n" +
                " * all tasks defined by a specific step in 'config.yaml'\n" +
                " * one or more steps\n" +
                "To start the complete pipeline as configured in 'config.yaml' " +
                "execute:\n" +
                "$ ./run-locally.py\n" +
                "To start a specific step execute:\n" +
                "$ ./run-locally.py <step_name>\n" +
                "To start a specific task execute:\n" +
                "$ ./run-locally.py <step_name/run_id>\n" +
                "The step_name is the name of an entry in the 'steps:' section " +
                "as defined in 'config.yaml'. A specific task is defined via " +
                "its task ID 'step_name/run_id'. A list of all task IDs is " + 
                "returned by running './status.py'.",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("--even-if-dirty",
                    dest="even_if_dirty",
                    action="store_true",
                    default=False,
                    help="Must be set if the local git repository " +
                    "contains uncommited changes. Otherwise the pipeline " +
                    "will not start.")

parser.add_argument("-s", "--step",
                    dest="step",
                    nargs='*',
                    default=list(),
                    type=str,
                    help="Can take multiple step names as input. A step name " +
                    "is the name of any entry in the 'steps:' section " +
                    "as defined in 'config.yaml'")

parser.add_argument("-t","--task",
                    dest="task",
                    nargs='*',
                    default=list(),
                    type=str,
                    help="Can take multiple task ID(s) as input. A task ID " +
                    "looks like ths 'step_name/run_id'. A list of all task IDs " +
                    "is returned by running './status.py'.")

args = parser.parse_args()

def main():
    p = pipeline.Pipeline(arguments=args)

    def handle_signal(signum, frame):
        print("Catching %s!" % process_pool.ProcessPool.SIGNAL_NAMES[signum])
        p.caught_signal = signum
        process_pool.ProcessPool.kill()
        
    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)

    task_list = copy.deepcopy(p.all_tasks_topologically_sorted)

    all_tasks = args.step + args.task

    if len(all_tasks) >= 1:
        # execute the specified tasks
        task_list = list()
        for task_id in sys.argv[1:]:
            task = p.task_for_task_id[task_id]
            task_list.append(task)
            
    # execute all tasks
    for task in task_list:
        basic_task_state = task.get_task_state_basic()
        if basic_task_state == p.states.FINISHED:
            sys.stderr.write("Skipping %s because it's already finished.\n" % task)
            continue
        if basic_task_state == p.states.READY:
            task.run()
        else:
            raise StandardError("Unexpected basic task state for %s: %s" % (task, basic_task_state))

if __name__ == '__main__':
    try:
        main()
    finally:
        # make sure all child processes get terminated
        process_pool.ProcessPool.kill_all_child_processes()
