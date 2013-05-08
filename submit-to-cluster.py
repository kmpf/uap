#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import copy
import os
import re
import subprocess
import yaml

def main():
    original_argv = copy.copy(sys.argv)

    p = pipeline.Pipeline()

    task_wish_list = None
    if len(sys.argv) > 1:
        task_wish_list = list()
        for _ in sys.argv[1:]:
            if '/' in _:
                task_wish_list.append(_)
            else:
                for task in p.all_tasks:
                    if str(task)[0:len(_)] == _:
                        task_wish_list.append(str(task))

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

    quotas = dict()
    quotas['default'] = 5

    # read quotas
    # -> for every step, a quota can be defined (with a default quota in place for steps
    #    which have no defined quota)
    # -> during submitting, there is a list of N previous job ids in which every item
    #    holds one of the previously submitted tasks
    if os.path.exists("quotas.yaml"):
        all_quotas = yaml.load(open("quotas.yaml", 'r'))
        hostname = subprocess.check_output(['hostname']).strip()
        for key in all_quotas.keys():
            if re.match(key, hostname):
                print("Applying quotas for " + hostname + ".")
                quotas = all_quotas[key]

    if not 'default' in quotas:
        raise StandardError("No default quota defined for this host.")

    quota_jids = {}
    quota_offset = {}

    def submit_task(task, dependent_tasks_in = []):
        dependent_tasks = copy.copy(dependent_tasks_in)

        step_name = next_task.step.get_step_name()
        if not step_name in quota_jids:
            size = quotas[step_name] if step_name in quotas else quotas['default']
            quota_jids[step_name] = [None for _ in range(size)]
            quota_offset[step_name] = 0

        quota_predecessor = quota_jids[step_name][quota_offset[step_name]]
        if quota_predecessor:
            dependent_tasks.append(quota_predecessor)

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
        args = ['./run-locally.py']
        if '--even-if-dirty' in original_argv:
            args.append('--even-if-dirty')
        args.append('"' + str(task) + '"')
        submit_script = submit_script.replace("#{COMMAND}", ' '.join(args))

        long_task_id_with_run_id = '%s_%s' % (str(task.step), task.run_id)
        long_task_id = str(task.step)
        short_task_id = long_task_id[0:15]

        qsub_args = ['qsub', '-N', short_task_id]
        qsub_args.append('-e')
        qsub_args.append(os.path.join(task.step.get_output_directory(), '.' + long_task_id_with_run_id + '.stderr'))
        qsub_args.append('-o')
        qsub_args.append(os.path.join(task.step.get_output_directory(), '.' + long_task_id_with_run_id + '.stdout'))

        # create the output directory if it doesn't exist yet
        # this is necessary here because otherwise, qsub will complain
        if not os.path.isdir(task.step.get_output_directory()):
            os.makedirs(task.step.get_output_directory())

        if len(dependent_tasks) > 0:
            qsub_args.append("-hold_jid")
            qsub_args.append(','.join(dependent_tasks))

        really_submit_this = True
        if task_wish_list:
            if not str(task) in task_wish_list:
                really_submit_this = False
        if really_submit_this:
            sys.stdout.write("Submitting task " + str(task) + " with " + str(task.step._cores) + " cores => ")
            process = subprocess.Popen(qsub_args, bufsize = -1, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
            process.stdin.write(submit_script)
            process.stdin.close()
            process.wait()
            response = process.stdout.read()
            job_id = re.search('Your job (\d+)', response).group(1)
            if job_id == None or len(job_id) == 0:
                raise StandardError("Error: We couldn't parse a job id from this:\n" + response)

            quota_jids[step_name][quota_offset[step_name]] = job_id
            quota_offset[step_name] = (quota_offset[step_name] + 1) % len(quota_jids[step_name])

            job_id_for_task[str(task)] = job_id

            print("%s (%s)" % (job_id, short_task_id))
            if len(dependent_tasks) > 0:
                print(" - with dependent tasks: " + ', '.join(dependent_tasks))
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

if __name__ == '__main__':
    main()
