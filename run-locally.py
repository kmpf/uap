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
