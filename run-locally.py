#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import pipeline
import unix_pipeline
import yaml

p = pipeline.Pipeline()

if len(sys.argv) > 1:
    # execute the specified task
    for task_id in sys.argv[1:]:
        task = p.task_for_task_id[task_id]
        task.run()
    exit(0)

task_list = copy.deepcopy(p.all_tasks)

try:
    while p.has_unfinished_tasks(task_list):
        task = p.pick_next_ready_task(task_list)
        task_list.remove(task)
        task.run()
finally:
    # make sure we kill all child processes on exit
    # TODO: Will this also kill child processes which have not been
    # forked from our script? If yes, that would be a problem.
    unix_pipeline.kill_all_child_processes()
