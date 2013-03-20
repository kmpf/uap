#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import pipeline
import yaml

p = pipeline.Pipeline()
task_list = copy.copy(p.all_tasks)

while p.has_unfinished_tasks(task_list):
    task = p.pick_next_ready_task(task_list)
    task_list.remove(task)
    task.run()
