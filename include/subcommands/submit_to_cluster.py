#!/usr/bin/env python

import sys
import datetime
import copy
import logging
import os
import re
import subprocess
import yaml

import abstract_step
import fscache
import pipeline
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

    if args.job_ids:
        ids = p.get_cluster_job_ids()
        print(' '.join(ids))
        return
    elif args.first_error:
        report_script = p.get_cluster_command('last_error')
        cmd = [os.path.join(args.uap_path, report_script)]
        if report_script == '':
            raise UAPError('A first error report script is not implemented '
                    'for the cluster "%s".' % p.get_cluster_type())
        elif not os.path.exists(cmd[0]):
            raise UAPError('The configured report script "%s" does not exist.'
                    % cmd[0])
        ids = p.get_cluster_job_ids()
        if len(ids)==0:
            raise UAPError('No jobs found.')
        cmd.extend(ids)
        try:
            process = subprocess.call(cmd)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise UAPError('The configured first error report script "%s" '
                        'closed with an error.' % report_script)
            else:
                raise e
        return

    task_wish_list = None
    if len(args.run) >= 1:
        task_wish_list = list()
        for _ in args.run:
            if '/' in _:
                task_wish_list.append(_)
            else:
                for task in p.all_tasks_topologically_sorted:
                    if str(task)[0:len(_)] == _:
                        task_wish_list.append(str(task))

    template_path = os.path.join(p.get_uap_path(),
                                 p.get_cluster_command('template'))
    try:
        template = open(template_path, 'r').read()
    except OSError:
        raise UAPError("Couldn't open %s." % template_path)

    steps_left = list()
    tasks_left = dict()
    for step_name in p.topological_step_order:
        tasks_left[step_name] = list()
        for task in p.tasks_in_step[step_name]:
            if task_wish_list is not None:
                if not str(task) in task_wish_list:
                    continue
            state = task.get_task_state()
            if state in [p.states.QUEUED, p.states.EXECUTING, p.states.FINISHED]:
                print("Skipping %s because it is already %s..." % (str(task), state.lower()))
                continue
            if step_name not in steps_left:
                steps_left.append(step_name)
            tasks_left[step_name].append(task)

    print("Now attempting to submit %d jobs..." % len(steps_left))

    quotas = dict()

    try:
        quotas['default'] = p.config['cluster']['default_job_quota']
        print("Set default quota to %s" % quotas['default'])
    except:
        print("No default quota defined in %s. Set default quota to '5'." %
              p.args.config.name)
        quotas['default'] = 5
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

        config_file_path = p.args.config.name
        command = [os.path.join(p.get_uap_path(), 'uap')]
        if p.args.debugging:
            command.append('--debugging')
        if p.args.verbose > 1:
            command.append('-' + 'v'*(p.args.verbose-1))
        command.extend([config_file_path, 'run-locally'])

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
        long_task_id_with_date = '%s_%%A_%%a_%s' % (step_name, now)

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

        submit_script_args.append(p.get_cluster_command('set_stderr'))
        submit_script_args.append(
            os.path.join(step.get_output_directory(),
                         '.' + long_task_id_with_date + '.stderr'))
        submit_script_args.append(p.get_cluster_command('set_stdout'))
        submit_script_args.append(
            os.path.join(step.get_output_directory(),
                         '.' + long_task_id_with_date + '.stdout'))

        if len(dependent_steps) > 0:
            submit_script_args += p.get_cluster_command_cli_option(
                'hold_jid',
                p.get_cluster_command('hold_jid_separator').join(dependent_steps))

        ##################
        # Submit the run #
        ##################
        print("Submitting step %s with %s cores per job"
              % (step_name, str(step._cores)))
        print("Submit command: %s" % " ".join(submit_script_args))
        print("=>")
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
        response = process.stdout.read()
        print("GOT A RESPONSE: %s" % response)
        job_id = re.search(
            p.get_cluster_command('parse_job_id'), response).group(1)

        if job_id == None or len(job_id) == 0:
            raise StandardError("Error: We couldn't parse a job_id from this:\n" + response)

        queued_ping_info = dict()
        queued_ping_info['step'] = step_name
        queued_ping_info['job_id'] = job_id
        queued_ping_info['submit_time'] = datetime.datetime.now()
        for task in tasks:
            queued_ping_info['run_id'] = task.run_id
            ping_file = step.get_run(task.run_id).get_queued_ping_file()
            if os.path.exists(ping_file+'.last'):
                os.unlink(ping_file+'.last')
            with open(ping_file, 'w') as f:
                f.write(yaml.dump(queued_ping_info, default_flow_style = False))

        print("%s (%s)" % (job_id, step_name))
        if len(dependent_steps) > 0:
            print(" - with dependent steps: " + ', '.join(dependent_steps))

        abstract_step.AbstractStep.fsc = fscache.FSCache()

    # After defining submit_task() let's walk through steps_left

    for step_name in steps_left:
        step = p.get_step(step_name)
        if not step_name in quotas.keys():
            quotas[step_name] = quotas['default']
            if step._options['_cluster_job_quota']:
                quotas[step_name] = step._options['_cluster_job_quota']
            print("Set job quota for %s to %s" % (step_name, quotas[step_name]))
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
                elif parent_state in [p.states.READY, p.states.WAITING]:
                    print("Cannot submit %s because a parent job "
                        "%s is %s when it should be queued, running, "
                        "or finished and the task selection as defined by "
                        "your command-line arguments do not request it to "
                        "submitted." % (step_name, parent_task, parent_state.lower()))
                    continue
        submit_step(step_name, parent_job_ids)
