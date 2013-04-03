#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import yaml

p = pipeline.Pipeline()

if len(sys.argv) > 1:
    # print a specific task
    task = p.task_for_task_id[sys.argv[1]]
    print(yaml.dump(task.step.get_run_info()[task.run_id], default_flow_style = False))
else:
    # print all tasks
    p.print_tasks()
