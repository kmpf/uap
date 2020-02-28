#!/usr/bin/env python
# encoding: utf-8

import sys
import os
from contextlib import closing
import logging
import pydoc
import string
from cStringIO  import StringIO
import yaml
from tqdm import tqdm
from subprocess import list2cmdline
import signal

import pipeline
from uaperrors import UAPError
import misc

'''
By default, this script displays information about all tasks of the pipeline
configured in 'config.yaml'. But the displayed information can be narrowed
down via command line options.

'''

logger = logging.getLogger("uap_logger")

def main(args):

    def exit(signum, frame):
        sys.exit(1)
    signal.signal(signal.SIGPIPE, exit)

    if args.sources:
        # print all sources (i. e. instances of AbstractSourceStep)
        args.no_tool_checks = True
        pipeline.Pipeline(arguments=args).print_source_runs()

    elif args.job_ids:
        args.no_tool_checks = True
        p = pipeline.Pipeline(arguments=args)
        ids = p.get_cluster_job_ids()
        print(' '.join(ids))
        return

    elif args.run and not args.details:
        # print run infos of one or more specific tasks
        p = pipeline.Pipeline(arguments=args)
        for task_id in args.run:
            parts = task_id.split('/')
            if len(parts) != 2:
                raise StandardError("Invalid task ID %s." % task_id)
            step_name = parts[0]
            run_id = parts[1]
            report = p.steps[step_name].get_run(run_id).as_dict()
            print(yaml.dump(report, Dumper=misc.UAPDumper,
                    default_flow_style = False))

    elif args.graph:
        p = pipeline.Pipeline(arguments=args)
        step_order = p.topological_step_order
        indents = [0 for _ in step_order]
        for index, line in enumerate(step_order):
            step_name = step_order[index]
            child_count = len(p.steps[step_name].children_step_names)
            indent = child_count * 2
            for _ in range(index + 1, len(step_order)):
                indents[_] += indent
            lines = list()
            for index, step_name in enumerate(step_order):
                lines.append(list(' ' * indents[index]))

        # draw horizontal line parts
        for index, step_name in enumerate(step_order):
            child_order = [_ for _ in step_order if _ in p.steps[step_name].children_step_names]
            for child_index, child in enumerate(child_order):
                x0 = indents[index] + 1 + child_index * 2
                y = step_order.index(child)
                x1 = indents[y]
                for x in range(x0, x1):
                    lines[y][x] = "─"
                lines[y][x0 - 1] = "└"

        # draw vertical line parts
        for index, step_name in enumerate(step_order):
            child_order = [_ for _ in step_order if _ in p.steps[step_name].children_step_names]
            for child_index, child in enumerate(child_order):
                x = indents[index] + child_index * 2
                y0 = index + 1
                y1 = step_order.index(child)
                for y in range(y0, y1):
                    lines[y][x] = "│"

        # write step names and states
        for index, line in enumerate(lines):
            original_step_name_label = ''
            step = p.steps[step_order[index]]
            if step.get_step_name() != step.get_step_type():
                original_step_name_label = ' (%s)' % step.get_step_type()

            if args.no_tool_checks:
                run_info = '%d runs' % len(step.get_runs())
            else:
                run_info = step.get_run_info_str(progress=True, do_hash=args.hash)
            line = "%s%s%s [%s]" % (''.join(line).replace("─└", "─┴"),
                    step.get_step_name(), original_step_name_label, run_info)
            print(line)

    elif args.details:
        p = pipeline.Pipeline(arguments=args)
        new_dest = p.config['destination_path']
        observed_states = set()
        tasks = list()
        for task_id in args.run:
            if task_id not in p.task_for_task_id.keys():
                raise UAPError('Task ID "%s" is not recognized.' % task_id)
            tasks.append(p.task_for_task_id[task_id])
        if not args.run:
            tasks = p.all_tasks_topologically_sorted
        for i, task in enumerate(tasks):
            state = task.get_task_state(do_hash=args.hash)
            sys.stdout.write('[%d/%d] ' % (i+1, len(tasks)))
            observed_states.add(state)
            if state == p.states.FINISHED:
                message = '%s is finished' % task
                if args.hash:
                    message += ' and sha256sum(s) correct'
                print(message)

            elif state == p.states.READY:
                print('%s has all inputs and is ready to start running' % task)

            elif state == p.states.EXECUTING:
                exec_ping_file = task.get_run().get_executing_ping_file()
                try:
                    with open(exec_ping_file, 'r') as buff:
                        info = yaml.load(buff,
                                Loader=yaml.FullLoader)
                except IOError as e:
                    if os.path.exists(exec_ping_file):
                        raise e
                    else:
                        print('%s is executing but seems to stop just now' % task)
                else:
                    print('%s is executing since %s' %
                            (task, info['start_time']))

            elif state == p.states.QUEUED:
                queued_ping_file = task.get_run().get_queued_ping_file()
                try:
                    with open(queued_ping_file, 'r') as buff:
                        info = yaml.load(buff,
                                Loader=yaml.FullLoader)
                except IOError as e:
                    if os.path.exists(queued_ping_file):
                        raise e
                    print('%s is queued but seems to stop just now' % task)
                else:
                    print('%s is queued with id %s since %s' %
                            (task, info['job_id'], info['submit_time']))

            elif state == p.states.VOLATILIZED:
                print('%s was volatilized and must be re-run if the data is '
                      'needed' % task)

            elif state == p.states.WAITING:
                parents = [str(parent) for parent in task.get_parent_tasks()
                        if parent.get_task_state() != p.states.FINISHED]
                print('%s is waiting for %s' % (task, parents))

            elif state == p.states.CHANGED:
                heading = 'changes in task %s' % task
                print(heading)
                print('-'*len(heading))
                run = task.get_run()
                anno_data = run.written_anno_data()
                has_date_change = False
                has_only_date_change = True
                if not anno_data:
                    has_only_date_change = False
                    print('No annotation file.')
                else:
                    changes = run.get_changes()
                    if changes:
                        has_only_date_change = False
                        print(yaml.dump(dict(changes), Dumper=misc.UAPDumper,
                                default_flow_style = False))
                    for line in run.file_changes(do_hash=args.hash,
                            report_correct=True):
                        if ' modification date after' in line:
                            has_date_change = True
                        else:
                            has_only_date_change = False
                        print(line)
                if has_date_change and has_only_date_change:
                    print("\nThis may be fixed with "
                         "'uap %s fix-problems --file-modification-date "
                         "--srsly'." %
                         args.config.name)
                print('')

            elif state == p.states.BAD:
                heading = 'errors of task %s' % task
                print(heading)
                print('-'*len(heading)+'\n')
                run = task.get_run()
                found_error = False
                anno_file = run.get_annotation_path()
                stale = run.is_stale()
                if stale:
                    print('Marked as executing but did not show activity '
                          'for %s.\n' % misc.duration_to_str(stale))
                anno_data = run.written_anno_data()
                if not anno_data:
                    print('The annotation file could not be read: %s.\n' %
                            anno_file)
                else:
                    failed = dict()
                    try:
                        procs = anno_data['pipeline_log'].get('processes', [])
                        for proc in procs:
                            if proc.get('exit_code', 0) == 0:
                                continue
                            failed[proc['name']] = {
                                'command':list2cmdline(proc['args']),
                                'exit code':proc['exit_code'],
                                'stderr':proc['stderr_copy']['tail']
                            }
                    except KeyError as e:
                        print('The annotation file "%s" seems badly '
                                'formated: %s\n' % (anno_file, e))
                    else:
                        host = anno_data.get('run', dict()).get('hostname', 'unknown')
                        time = anno_data.get('end_time', 'unknown')
                        print('host: %s' % host)
                        print('time: %s' % time)
                        print('')
                    if failed:
                        found_error = True
                        print('## FAILED COMMANDS ##')
                        print(yaml.dump(failed, Dumper=misc.UAPDumper,
                                default_flow_style = False))
                    else:
                        print('No failed commands found in the annotation file.\n')
                    run_data = anno_data.get('run', [])
                    if 'error' in run_data:
                        found_error = True
                        print('## ERROR ##\n%s\n' %  run_data['error'])
                    if not found_error:
                        print('No errors found.')
                        print("Run 'uap %s fix-problems --first-error' to investigate.'"
                                % args.config.name)
                        print('')

            else:
                print('%s has unknown state "%s"' % (task, state.lower()))

        if p.states.CHANGED in observed_states:
            print("If you want to force overwrite of the changed tasks, run\n"
                  "'uap %s run-locally --force' or 'uap %s submit-to-cluster --force'." %
                  (args.config.name, args.config.name))
            print("If you want to ignore the changes and consider the tasks finished, run\n"
                  "'uap %s run-locally --irgnore' or 'uap %s submit-to-cluster --irgnore'." %
                  (args.config.name, args.config.name))

        # now check ping files and print some warnings and instructions if
        # something's fishy
        p.check_ping_files(print_more_warnings = True if args.verbose > 0 else False)

        # Now check whether we can volatilize files, but don't do it.
        p.check_volatile_files()

    else:
        # print all runs
        '''
        prints a summary of all runs, indicating whether each run is
        - ``[r]eady``
        - ``[w]aiting``
        - ``[q]ueued``
        - ``[e]xecuting``
        - ``[b]ad``
        - ``[f]inished``
        - ``[c]hanged``
        - ``[v]olatilized``
        '''
        p = pipeline.Pipeline(arguments=args)
        output = list()
        tasks_for_status = {}
        tasks = p.all_tasks_topologically_sorted

        task_iter = tqdm(tasks, desc='tasks')
        try:
            for task in task_iter:
                state = task.get_task_state(do_hash=args.hash)
                if not state in tasks_for_status:
                    tasks_for_status[state] = list()
                tasks_for_status[state].append(task)
        except:
            task_iter.close()
            raise

        for status in p.states.order:
            if not status in tasks_for_status:
                continue
            heading = "%s runs" % string.capwords(status)
            output.append(heading)
            output.append('-' * len(heading))
            if args.summarize:
                step_count = dict()
                step_order = list()
                for task in tasks_for_status[status]:
                    if not str(task.step) in step_count:
                        step_count[str(task.step)] = 0
                        step_order.append(str(task.step))
                    step_count[str(task.step)] += 1
                for step_name in step_order:
                    output.append("[%s]%4d %s"
                                 % (status.lower()[0],
                                    step_count[step_name],
                                    step_name)
                             )
            else:
                for task in tasks_for_status[status]:
                    output.append("[%s] %s"
                                 % (
                                     status[0].lower(),
                                     task))
                output.append('')
        output.append("runs: %d total, %s"
                      % (len(p.all_tasks_topologically_sorted),
                         ', '.join(["%d %s" % (
                             len(tasks_for_status[_]),
                             _.lower()) for _ in p.states.order \
                                    if _ in tasks_for_status])))
        pydoc.pager("\n".join(output))

        if p.states.CHANGED in tasks_for_status.keys() \
        or p.states.BAD in tasks_for_status.keys() \
        or p.states.WAITING in tasks_for_status.keys():
            print("\nRun 'uap %s status --details' to inspect states." %
                    args.config.name)

        # now check ping files and print some warnings and instructions if
        # something's fishy
        p.check_ping_files(print_more_warnings = True if args.verbose > 0 else False)

        # Now check whether we can volatilize files, but don't do it.
        p.check_volatile_files()

