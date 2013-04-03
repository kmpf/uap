#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import pipeline
import yaml

p = pipeline.Pipeline()

if '--run-this' in sys.argv:
    # execute the specified task
    task_id = sys.argv[sys.argv.index('--run-this') + 1]
    task = p.task_for_task_id[task_id]
    task.run()
    exit(0)

task_list = copy.deepcopy(p.all_tasks)

while p.has_unfinished_tasks(task_list):
    task = p.pick_next_ready_task(task_list)
    task_list.remove(task)
    task.run()
