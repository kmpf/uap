#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import yaml

p = pipeline.Pipeline()

tasks_left = []

# a hash of files which are already there or will be there once submitted jobs
# have completed
file_hash = {}

for task in p.all_tasks:
    state = task.get_task_state()
    if state != p.states.FINISHED:
        tasks_left.append(task)
    else:
        for path in task.output_files():
            if not path in file_hash:
                file_hash[path] = []
            file_hash[path].append(str(task))

job_id_for_task = {}

def submit_task(task, dependent_tasks = None):
    print("Submitting task " + str(task) + " with " + str(task.step._cores) + " cores.")
    if dependent_tasks != None:
        print("...with dependent tasks: " + str(dependent_tasks))
    for path in task.output_files():
        if not path in file_hash:
            file_hash[path] = []
        file_hash[path].append(str(task))

    tasks_left.remove(task)
    job_id_for_task[str(task)] = 'JID_' + str(len(job_id_for_task) + 1)

# first submit all tasks which are ready as per the file system
while len(tasks_left) > 0:
    next_task = None
    for task in tasks_left:
        if task.get_task_state() == p.states.READY:
            next_task = task
            break

    if next_task == None:
        break

    submit_task(next_task)

while len(tasks_left) > 0:
    next_task = None
    dependent_tasks = set()

    for task in tasks_left:
        # see which input files are required by the task and whether they are already in file_hash
        can_start = True
        dependent_tasks = set()
        for path in task.input_files():
            if path in file_hash:
                for t in file_hash[path]:
                    dependent_tasks.add(t)
            else:
                can_start = False
                break
        if can_start:
            next_task = task
            break

    if next_task == None:
        print("Error: Unable to find next task while there are still tasks left, giving up :(")
        exit(1)

    if len(dependent_tasks) == 0:
        dependent_tasks = None
    else:
        dependent_tasks = [job_id_for_task[_] for _ in dependent_tasks]

    submit_task(next_task, dependent_tasks)

print("All tasks submitted successfully.")
