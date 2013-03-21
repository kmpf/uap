#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import pipeline
import yaml

run_mode = pipeline.Pipeline.run_modes.FULL
if sys.argv[1] == '--dry-run':
    run_mode = pipeline.Pipeline.run_modes.DRY_RUN
elif sys.argv[1] == '--test-run':
    run_mode = pipeline.Pipeline.run_modes.TEST_RUN

p = pipeline.Pipeline(run_mode)

task_list = copy.copy(p.all_tasks)

while p.has_unfinished_tasks(task_list):
    task = p.pick_next_ready_task(task_list)
    task_list.remove(task)
    task.run()
