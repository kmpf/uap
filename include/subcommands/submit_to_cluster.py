#!/usr/bin/env python

import sys
import datetime
import copy
import logging
import os
import re
import subprocess
import yaml
from tqdm import tqdm

import abstract_step
import fscache
import pipeline
import submit_to_cluster_legacy
from uaperrors import UAPError
'''
By default, this script submits all tasks to a compute cluster via a
submit script. The list of tasks can be narrowed down by specifying a step name
(in which case all runs of this steps will be considered) or individual
tasks (step_name/run_id).

At this point, we have a task wish list.

This task wish list is now processed one by one (in topological order):

- determine status of current task
- if it's finished or running or queued, skip it
- if it's ready, submit it
- if it's waiting, there is at least one parent task which is not yet finished
- determine all parent tasks which are not finished yet, for each parent:
  - if it's running or queued, determine its job_id from the queued-ping file
  - if it's ready or waiting, we must skip this task because it should have been
    enqueued a couple of iterations ago (this is because a selection has been
    made and there are unfinished, non-running, unqueued dependencies)
  - now add all these collected job_ids to the submission via -hold_jid
    (or whatever the argument is for the cluster used)
'''

logger = logging.getLogger("uap_logger")

def main(args):
    if args.legacy is True:
        submit_to_cluster_legacy.main(args)
        return
    p = pipeline.Pipeline(arguments=args)

    template_path = os.path.join(p.get_uap_path(),
                                 p.get_cluster_command('template'))
    try:
        template = open(template_path, 'r').read()
    except OSError:
        raise UAPError("Couldn't open %s." % template_path)

    steps_left = list()
    tasks_left = dict()
    skip_message = list()
    skipped_tasks = dict()
    iter_steps = tqdm(p.topological_step_order, desc='step states')
    try:
        for step_name in iter_steps:
            tasks_left[step_name] = list()
            for task in p.tasks_in_step[step_name]:
                if p.task_wish_list:
                    if not str(task) in p.task_wish_list:
                        continue
                state = task.get_task_state()
                if state in [p.states.QUEUED, p.states.EXECUTING, p.states.FINISHED]:
                    skipped_tasks.setdefault(state, set())
                    skipped_tasks[state].add(task)
                    continue
                if state == p.states.VOLATILIZED and not p.task_wish_list:
                    skipped_tasks.setdefault(state, set())
                    skipped_tasks[state].add(task)
                    continue
                if state == p.states.CHANGED and args.ignore:
                    skipped_tasks.setdefault(state, set())
                    skipped_tasks[state].add(task)
                    skip_message.append("Skipping %s because it's changes are "
                        "ignored." % task)
                    continue
                if state == p.states.CHANGED and not args.force:
                    raise UAPError("Task %s has changed. "
                            "Run 'uap %s status --details' to see what changed or "
                            "'uap %s submit-to-cluster --force' to force overwrite "
                            "of the results." %
                            (task, args.config.name, args.config.name))
                tasks_left[step_name].append(task)
                if step_name not in steps_left:
                    steps_left.append(step_name)
    except:
        iter_steps.close()
        raise

    for state, tasks in skipped_tasks.items():
        skip_message.append("Skipping %d task(s) because they are %s..." %
                (len(tasks), state.lower()))
    if p.states.VOLATILIZED in skipped_tasks.keys():
        skip_message.append("Request reproduction of volitilized tasks by "
                            "requesting it explicitly, e.g., with 'uap %s "
                            "submit-to-cluster <step name>'." %
                            p.args.config.name)

    quotas = dict()
    for line in skip_message:
        print(line)

    try:
        quotas['default'] = p.config['cluster']['default_job_quota']
    except:
        quotas['default'] = 0
    # read quotas
    # -> for every step, a quota can be defined (with a default quota
    #    in place for steps which have no defined quota)
    # -> during submission, there is a list of N previous job ids in
    #    which every item holds one of the previously submitted tasks

    def submit_step(step_name, dependent_steps = []):
        '''
        This method reads and modifies the necessary submit script for a given
        step. It applies job quotas. Finally, it starts the submit command.
        '''
        step = p.get_step(step_name)

        ##########################
        # Assemble submit script #
        ##########################

        submit_script = copy.copy(template)
        # Here we need to replace the place holders with the stuff defined
        # in our config file
        #

        # Get the values for the placeholders:
        # SUBMIT_OPTIONS
        placeholder_values = dict()
        if step._options['_cluster_submit_options']:
            placeholder_values['#{SUBMIT_OPTIONS}'] = step._options[
                '_cluster_submit_options']
        else:
            placeholder_values['#{SUBMIT_OPTIONS}'] = p.config['cluster']\
                                                      ['default_submit_options']
        # PRE_JOB_COMMAND
        if step._options['_cluster_pre_job_command']:
            placeholder_values['#{PRE_JOB_COMMAND}'] = step._options[
            '_cluster_pre_job_command']
        else:
            placeholder_values['#{PRE_JOB_COMMAND}'] = p.config['cluster']\
                                                       ['default_pre_job_command']
        # POST_JOB_COMMAND
        if step._options['_cluster_post_job_command']:
            placeholder_values['#{POST_JOB_COMMAND}'] = step._options[
                '_cluster_post_job_command']
        else:
            placeholder_values['#{POST_JOB_COMMAND}'] = p.config['cluster']\
                                                        ['default_post_job_command']

        # Replace placeholders with their values
        for placeholder, value in placeholder_values.items():
            submit_script = submit_script.replace(placeholder, value)

        task_names = [str(task) for task in tasks_left[step_name]]
        submit_script = submit_script.replace("#{ARRAY_JOBS}",
                " ".join("'" + task + "'" for task in task_names))
        submit_script = submit_script.replace("#{CORES}", str(step._cores))
        submit_script = submit_script.replace("#{UAP_CONFIG}", yaml.dump(p.config))

        command = ['exec', os.path.join(p.get_uap_path(), 'uap'), '-vv']
        if p.args.debugging:
            command.append('--debugging')
        command.extend(['<(cat <&123)', 'run-locally'])
        if p.args.force:
            command.append('--force')

        task_id = p.get_cluster_command('array_task_id')
        command.append('"${array_jobs[$' + task_id + ']}"')

        submit_script = submit_script.replace("#{COMMAND}", ' '.join(command))

        # create the output directory if it doesn't exist yet
        tasks = tasks_left[step_name]
        for task in tasks:
            run_output_dir = task.get_run().get_output_directory()
            if not os.path.isdir(run_output_dir):
                os.makedirs(run_output_dir)

        ###########################
        # Assemble submit command #
        ###########################
        now = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        aoi = p.get_cluster_command('array_out_index')
        if aoi and aoi != '':
            long_task_id_with_date = '_'.join([step_name, aoi, now])
        else:
            long_task_id_with_date = '_'.join([step_name, now])

        submit_script_args = [p.get_cluster_command('submit')]
        size = len(tasks)
        if quotas[step_name] == 0:
            submit_script_args += p.get_cluster_command_cli_option('array_job',
                    str(size-1))
        else:
            submit_script_args += p.get_cluster_command_cli_option(
                    'array_job_wquota',
                    (str(size-1), str(quotas[step_name])))
        submit_script_args += p.get_cluster_command_cli_option(
            'set_job_name', step_name)

        out_file = os.path.join(step.get_output_directory(),
                '.' + long_task_id_with_date + '.uapout')
        submit_script_args.append(p.get_cluster_command('set_stderr'))
        submit_script_args.append(out_file)
        submit_script_args.append(p.get_cluster_command('set_stdout'))
        submit_script_args.append(out_file)

        if len(dependent_steps) > 0:
            submit_script_args += p.get_cluster_command_cli_option(
                'hold_jid',
                p.get_cluster_command('hold_jid_separator').join(dependent_steps))

        ##################
        # Submit the run #
        ##################
        sys.stdout.write(" with %s cores per job" % str(step._cores))
        if quotas[step_name] != 0:
            sys.stdout.write(", quota %s" % quotas[step_name])
        else:
            sys.stdout.write(", no quota")
        if dependent_steps:
            sys.stdout.write(", dependencies %s" % ', '.join(dependent_steps))
        else:
            sys.stdout.write(", no dependencies")
        # Store submit script in the run_output_dir
        submit_script_path = step.get_submit_script_file()
        with open(submit_script_path, 'w') as f:
            f.write(submit_script)

        process = None
        try:
            process = subprocess.Popen(
                submit_script_args,
                bufsize = -1,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise StandardError("Unable to launch %s. Maybe " %
                                    p.get_cluster_command('submit') +
                                    "you are not executing this script on " +
                                    "the cluster")
            else:
                raise e
        process.stdin.write(submit_script)
        process.stdin.close()
        process.wait()
        response = process.stdout.read().strip()
        job_id = re.search(
            p.get_cluster_command('parse_job_id'), response)
        if not job_id:
            raise UAPError('Got unexpected from %s resposne: %s' %
                    (p.get_cluster_type(), response))
        else:
            job_id = job_id.group(1)
        sys.stdout.write(" and job id %s.\n" % job_id)

        if job_id == None or len(job_id) == 0:
            raise StandardError("Error: We couldn't parse a job_id from this:\n" + response)

        queued_ping_info = dict()
        queued_ping_info['step'] = step_name
        queued_ping_info['job_id'] = job_id
        queued_ping_info['submit_time'] = datetime.datetime.now()
        for task in tasks:
            queued_ping_info['run_id'] = task.run_id
            ping_file = step.get_run(task.run_id).get_queued_ping_file()
            if os.path.exists(ping_file+'.bad'):
                os.unlink(ping_file+'.bad')
            with open(ping_file, 'w') as f:
                f.write(yaml.dump(queued_ping_info, default_flow_style = False))

        task.get_run().reset_fsc()

    # After defining submit_task() let's walk through steps_left

    for step_num, step_name in enumerate(steps_left):
        step = p.get_step(step_name)
        if not step_name in quotas.keys():
            quotas[step_name] = quotas['default']
            if step._options['_cluster_job_quota']:
                quotas[step_name] = step._options['_cluster_job_quota']
        parents = step.get_dependencies()
        parent_job_ids = set()
        for parent in parents:
            for parent_task in p.tasks_in_step[parent.get_step_name()]:
                parent_state = parent_task.get_task_state()
                if parent_state in [p.states.EXECUTING, p.states.QUEUED]:
                    # determine job_id from YAML queued ping file
                    parent_job_id = None
                    parent_queued_ping_path = parent.get_run(parent_task.run_id).get_queued_ping_file()
                    try:
                        parent_info = yaml.load(open(parent_queued_ping_path), Loader=yaml.FullLoader)
                        parent_job_ids.add(parent_info['job_id'])
                    except:
                        print("Couldn't determine job_id of %s while trying to load %s." %
                            (parent_task, parent_queued_ping_path))
                        raise
                elif parent_state in [p.states.READY, p.states.WAITING, p.states.BAD, p.states.CHANGED]:
                    print("Cannot submit %s because its "
                        "parent %s is %s when it should be queued, running, "
                        "or finished." % (task, parent_task, parent_state.lower()))
                    continue
        sys.stdout.write("[%d/%d][%s] %s job" %
                (step_num+1, len(steps_left), step_name, p.get_cluster_type()))
        submit_step(step_name, parent_job_ids)
        step.reset_run_caches()
        sys.stdout.flush()
