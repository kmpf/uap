#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import pipeline
import process_pool
import yaml

def main():
    p = pipeline.Pipeline()

    if len(sys.argv) > 1:
        # execute the specified tasks
        for task_id in sys.argv[1:]:
            task = p.task_for_task_id[task_id]
            task.run()
    else:
        # execute all tasks
        task_list = copy.deepcopy(p.all_tasks)
        while p.has_unfinished_tasks(task_list):
            task = p.pick_next_ready_task(task_list)
            task_list.remove(task)
            task.run()

if __name__ == '__main__':
    try:
        main()
    finally:
        # make sure all child processes get terminated
        process_pool.ProcessPool.kill_all_child_processes()