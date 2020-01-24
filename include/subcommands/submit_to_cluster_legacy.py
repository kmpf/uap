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
    p = pipeline.Pipeline(arguments=args)

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

    tasks_left = []

    template_path = os.path.join(p.get_uap_path(),
                                 p.get_cluster_command('template'))
    try:
        template = open(template_path, 'r').read()
    except OSError:
        logger.error("Couldn't open %s." % template_path)
        sys.exit(1)

    for task in p.all_tasks_topologically_sorted:
        if task_wish_list is not None:
            if not str(task) in task_wish_list:
                continue
        tasks_left.append(task)

    print("Now attempting to submit %d jobs..." % len(tasks_left))

    quotas = dict()

    try:
        quotas['default'] = p.config['cluster']['default_job_quota']
        print("Set default quota to %s" % quotas['default'])
    except:
        print("No default quota defined in %s. Set default quota to '5'." %
              p.get_config_filepath())
        quotas['default'] = 5
    # read quotas
    # -> for every step, a quota can be defined (with a default quota
    #    in place for steps which have no defined quota)
    # -> during submission, there is a list of N previous job ids in
    #    which every item holds one of the previously submitted tasks

    quota_jids = {}
    quota_offset = {}

    def submit_task(task, dependent_tasks_in = []):
        '''
        This method reads and modifies the necessary submit script for a given
        task. It applies job quotas. Finally, it starts the submit command.
        '''
        dependent_tasks = copy.copy(dependent_tasks_in)

        step = task.get_step()
        step_name = task.step.get_step_name()
        size = quotas[step_name] if step_name in quotas else quotas['default']
        if not step_name in quota_jids:
            quota_jids[step_name] = [None for _ in range(size)]
            quota_offset[step_name] = 0

        if size > 0:
            quota_predecessor = quota_jids[step_name][quota_offset[step_name]]
            if quota_predecessor:
                dependent_tasks.append(quota_predecessor)

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

        submit_script = submit_script.replace("#{ARRAY_JOBS}", '')
        submit_script = submit_script.replace("#{CORES}", str(task.step._cores))

        # todo: set email address

        config_file_path = args.config.name
        command = ['uap', config_file_path, 'run-locally']

        command.append('"' + str(task) + '"')

        submit_script = submit_script.replace("#{COMMAND}", ' '.join(command))

        # create the output directory if it doesn't exist yet
        # this is necessary here because otherwise, qsub will complain
        run_output_dir = task.get_run().get_output_directory()
        if not os.path.isdir(run_output_dir):
            os.makedirs(run_output_dir)

        ###########################
        # Assemble submit command #
        ###########################
        now = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        long_task_id_with_run_id = '%s_%s_%s' % (str(task.step),
                                                 task.run_id, now)
        long_task_id = str(task.step)

        submit_script_args = [p.get_cluster_command('submit')]
        submit_script_args += p.get_cluster_command_cli_option(
            'set_job_name', long_task_id)

        submit_script_args.append(p.get_cluster_command('set_stderr'))
        submit_script_args.append(
            os.path.join(task.get_run().get_output_directory(),
                         '.' + long_task_id_with_run_id + '.stderr'))
        submit_script_args.append(p.get_cluster_command('set_stdout'))
        submit_script_args.append(
            os.path.join(task.get_run().get_output_directory(),
                         '.' + long_task_id_with_run_id + '.stdout'))

        if len(dependent_tasks) > 0:
            submit_script_args += p.get_cluster_command_cli_option(
                'hold_jid',
                p.get_cluster_command('hold_jid_separator').join(dependent_tasks))

        ########################################
        # Submit the run if everything is fine #
        ########################################
        really_submit_this = True
        if task_wish_list:
            if not str(task) in task_wish_list:
                really_submit_this = False
        if task.get_task_state() == p.states.EXECUTING:
            print("Skipping %s because it is already executing." % str(task))
            really_submit_this = False
        if really_submit_this:
            print("Submitting task %s with %s cores"
                  % (str(task), str(task.step._cores)))
            print("Submit command: %s" % " ".join(submit_script_args))
            print("=>")
            # Store submit script in the run_output_dir
            submit_script_path = task.get_run().get_submit_script_file()
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
            queued_ping_info['step'] = str(task.step)
            queued_ping_info['run_id'] = task.run_id
            queued_ping_info['job_id'] = job_id
            queued_ping_info['submit_time'] = datetime.datetime.now()
            with open(task.get_step().get_run(task.run_id).get_queued_ping_file(), 'w') as f:
                f.write(yaml.dump(queued_ping_info, default_flow_style = False))

            if size > 0:
                quota_jids[step_name][quota_offset[step_name]] = job_id
                quota_offset[step_name] = (quota_offset[step_name] + 1) % len(quota_jids[step_name])

            print("%s (%s)" % (job_id, long_task_id))
            if len(dependent_tasks) > 0:
                print(" - with dependent tasks: " + ', '.join(dependent_tasks))

            abstract_step.AbstractStep.fsc = fscache.FSCache()

    # After defining submit_task() let's walk through tasks_left

    for task in tasks_left:
        # Update quotas dict for current step if necessary
        step_name = task.step.get_step_name()
        if step_name not in quotas.keys():
            step = task.get_step()
            quotas[step_name] = quotas['default']
            if step._options['_cluster_job_quota']:
                quotas[step_name] = step._options['_cluster_job_quota']

            print("Set job quota for %s to %s" % (step_name, quotas[step_name]))
        state = task.get_task_state()
        if state in [p.states.QUEUED, p.states.EXECUTING, p.states.FINISHED]:
            print("Skipping %s because it is already %s..." % (str(task), state.lower()))
            continue
        if state == p.states.READY:
            submit_task(task)
        if state == p.states.WAITING:
            skip_this = False
            parent_job_ids = list()
            for parent_task in task.get_parent_tasks():
                parent_state = parent_task.get_task_state()
                if parent_state in [p.states.EXECUTING, p.states.QUEUED]:
                    # determine job_id from YAML queued ping file
                    parent_job_id = None
                    parent_queued_ping_path = parent_task.get_step().get_run(parent_task.run_id).get_queued_ping_file()
                    try:
                        parent_info = yaml.load(open(parent_queued_ping_path), Loader=yaml.FullLoader)
                        parent_job_ids.append(parent_info['job_id'])
                    except:
                        print("Couldn't determine job_id of %s while trying to load %s." %
                            (parent_task, parent_queued_ping_path))
                        raise
                elif parent_state in [p.states.READY, p.states.WAITING]:
                    skip_this = True
                    print("Cannot submit %s because its "
                        "parent %s is %s when it should be queued, running, "
                        "or finished and the task selection as defined by "
                        "your command-line arguments do not request it to "
                        "submitted." % (task, parent_task, parent_state.lower()))
            if not skip_this:
                submit_task(task, parent_job_ids)
