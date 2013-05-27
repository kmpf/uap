#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import misc
import os
import pipeline
import process_pool
import signal
import socket
import yaml

def main():
    p = pipeline.Pipeline()
    
    def handle_signal(signum, frame):
        print("Catching %s!" % process_pool.ProcessPool.SIGNAL_NAMES[signum])
        p.caught_signal = signum
        process_pool.ProcessPool.kill()
        
    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)

    task_list = copy.deepcopy(p.all_tasks_topologically_sorted)

    if len(sys.argv) > 1:
        # execute the specified tasks
        task_list = list()
        for task_id in sys.argv[1:]:
            task = p.task_for_task_id[task_id]
            task_list.append(task)
            
    # execute all tasks
    for task in task_list:
        if task.get_task_state() == p.states.FINISHED:
            continue
        if task.get_task_state() in [p.states.READY, p.states.QUEUED]:
            task.run()
        else:
            raise StandardError("Unexpected task state for %s: %s" % (task, task.get_task_state()))

if __name__ == '__main__':
    try:
        main()
    finally:
        # make sure all child processes get terminated
        process_pool.ProcessPool.kill_all_child_processes()
