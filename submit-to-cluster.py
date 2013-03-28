#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import copy
import re
import subprocess
import yaml

original_argv = copy.copy(sys.argv)

p = pipeline.Pipeline()

if '--run-this' in sys.argv:
    # execute the specified task
    task_id = sys.argv[sys.argv.index('--run-this') + 1]
    task = p.task_for_task_id[task_id]
    task.run()
    exit(0)

tasks_left = []

# a hash of files which are already there or will be there once submitted jobs
# have completed
file_hash = {}

template = ''
with open('qsub-template.sh', 'r') as f:
    template = f.read()

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

def submit_task(task, dependent_tasks = []):
    print("Submitting task " + str(task) + " with " + str(task.step._cores) + " cores.")
    if len(dependent_tasks) > 0:
        print("...with dependent tasks: " + str(dependent_tasks))
    for path in task.output_files():
        if not path in file_hash:
            file_hash[path] = []
        file_hash[path].append(str(task))

    submit_script = copy.copy(template)
    submit_script = submit_script.replace("#{CORES}", str(task.step._cores))
    email = 'nobody@example.com'
    if 'email' in p.config:
        email = p.config['email']
    submit_script = submit_script.replace("#{EMAIL}", email)
    args = copy.copy(original_argv)
    args.append("--run-this")
    args.append('"' + str(task) + '"')
    submit_script = submit_script.replace("#{COMMAND}", ' '.join(args))

    temp = str(task).split('/')
    for _ in range(len(temp) - 1):
        temp[_] = temp[_][0]
    short_task_id = (''.join(temp[0:-1]) + '_' + temp[-1])[0:15]

    qsub_args = ['qsub', '-N', short_task_id]
    qsub_args.append('-e')
    qsub_args.append(os.path.join(task.step.get_output_directory(), '.' + short_task_id + '.stderr'))
    qsub_args.append('-o')
    qsub_args.append(os.path.join(task.step.get_output_directory(), '.' + short_task_id + '.stdout'))
    if len(dependent_tasks) > 0:
        qsub_args.append("-hold_jid")
        qsub_args.append(','.join(dependent_tasks))

    process = subprocess.Popen(qsub_args, bufsize = -1, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    process.stdin.write(submit_script)
    process.stdin.close()
    process.wait()
    response = process.stdout.read()
    job_id = re.search('Your job (\d+)', response).group(1)
    if job_id == None or len(job_id) == 0:
        raise StandardError("Error: We couldn't parse a job id from this:\n" + response)

    job_id_for_task[str(task)] = job_id
    tasks_left.remove(task)

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

    jids = []
    if len(dependent_tasks) > 0:
        for _ in dependent_tasks:
            if _ in job_id_for_task:
                jids.append(job_id_for_task[_])

    submit_task(next_task, jids)

print("All tasks submitted successfully.")
