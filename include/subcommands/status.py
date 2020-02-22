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

    elif args.run and not args.hash:
        # print run infos of one or more specific tasks
        args.no_tool_checks = True
        p = pipeline.Pipeline(arguments=args)
        for task_id in args.run:
            parts = task_id.split('/')
            if len(parts) != 2:
                raise StandardError("Invalid task ID %s." % task_id)
            step_name = parts[0]
            run_id = parts[1]
            report = p.steps[step_name].get_run(run_id).as_dict()
            print(yaml.dump(report, default_flow_style = False))

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

        for index, _ in enumerate(lines):
            original_step_name_label = ''
            step = p.steps[step_order[index]]
            if step.get_step_name() != step.get_step_type():
                original_step_name_label = ' (%s)' % step.get_step_type()

#            sys.stderr.write("step_name: %s\n" % step.get_step_name())
#            sys.stderr.write("index: %d\n" % index)
#            sys.stderr.write("still alive!\n")
#            sys.stderr.write("run_info: %s\n" % step.get_run_info_str())

            if args.no_tool_checks:
                line = "%s%s%s" % (''.join(_).replace("─└", "─┴"),
                        step.get_step_name(), original_step_name_label)
            else:
                line = "%s%s%s [%s]" % (''.join(_).replace("─└", "─┴"),
                        step.get_step_name(), original_step_name_label,
                        step.get_run_info_str(progress=True))
            print(line)
        # now check ping files and print some warnings and instructions if
        # something's fishy
        p.check_ping_files(print_more_warnings = True if args.verbose > 0 else False)

        # Now check whether we can volatilize files, but don't do it.
        p.check_volatile_files()
    elif args.hash:
        p = pipeline.Pipeline(arguments=args)
        new_dest = p.config['destination_path']
        tasks = list()
        for task_id in args.run:
            if task_id not in p.task_for_task_id.keys():
                raise UAPError('Task ID "%s" is not recognized.' % task_id)
            tasks.append(p.task_for_task_id[task_id])
        if not args.run:
            tasks = p.all_tasks_topologically_sorted
        for task in tasks:
            state = task.get_task_state()
            if state in [p.states.READY, p.states.EXECUTING, p.states.BAD]:
                print('%s is %s' % (task, state.lower()))
            elif state == p.states.WAITING:
                parents = [str(parent) for parent in task.get_parent_tasks()
                        if parent.get_task_state() != p.states.FINISHED]
                print('%s is %s for %s' % (task, state.lower(), parents))
            elif state in [p.states.FINISHED, p.states.CHANGED]:
                title = '%s is %s and ' % (task, state.lower())
                sys.stdout.write(title)
                header = 'has changed files'
                header += '\n' + '-'*len(title+header) + '\n'
                good = True
                for bad_file in task.get_run().file_changes(do_hash=True):
                    if bad_file is None or isinstance(bad_file, bool):
                        break
                    if good:
                        sys.stdout.write(header)
                        good = False
                    print(bad_file)
                if bad_file is None:
                    print('has no output files')
                elif not bad_file:
                    print('the sha256sum(s) correct')
                else:
                    print('')
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
        '''
        p = pipeline.Pipeline(arguments=args)
        output = list()
        tasks_for_status = {}
        tasks = p.all_tasks_topologically_sorted

        for task in tqdm(tasks, desc='tasks'):
            state = task.get_task_state()
            if not state in tasks_for_status:
                tasks_for_status[state] = list()
            tasks_for_status[state].append(task)

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

        if p.states.WAITING in tasks_for_status.keys():
            if args.details:
                print('')
                for task in tasks_for_status[p.states.WAITING]:
                    parents = [str(parent) for parent in task.get_parent_tasks()
                            if parent.get_task_state() != p.states.FINISHED]
                    print('%s is waiting for %s' % (task, parents))

        if p.states.CHANGED in tasks_for_status.keys():
            if args.details:
                print('')
                for task in tasks_for_status[p.states.CHANGED]:
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
                            print(yaml.dump(dict(changes)))
                        for line in run.file_changes():
                            if line is None or isinstance(line, bool):
                                break
                            if ' was changed after' in line:
                                has_date_change = True
                            else:
                                has_only_date_change = False
                            print(line)
                    if has_date_change and has_only_date_change:
                        print("\nThis may be fixed with "
                             "'uap %s fix-problems --file-modification-date'." %
                             p.args.config.name)
                    print('')
            else:
                print("Some tasks changed. Run 'uap %s status --details' to see the details." %
                        p.args.config.name)
            print("If you want to force overwrite of the changed tasks, run\n"
                  "'uap %s run-locally --force' or 'uap %s submit-to-cluster --force'." %
                  (p.args.config.name, p.args.config.name))
            print("If you want to ignore the changes and consider the tasks finished, run\n"
                  "'uap %s run-locally --irgnore' or 'uap %s submit-to-cluster --irgnore'." %
                  (p.args.config.name, p.args.config.name))

        if p.states.BAD in tasks_for_status.keys():
            if args.details:
                print('')
                for task in tasks_for_status[p.states.BAD]:
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
                            procs = anno_data['pipeline_log']['processes']
                            for proc in procs:
                                if proc.get('exit_code', 0) == 0:
                                    continue
                                failed[proc['name']] = {
                                    'command':list2cmdline(proc['args']),
                                    'exit code':proc['exit_code'],
                                    'stderr':proc['stderr_copy']['tail']
                                }
                            run_data = anno_data['run']
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
                            print('#### failed commands')
                            print(yaml.dump(failed))
                        else:
                            print('No failed commands found in the annotation file.\n')
                        if 'error' in run_data:
                            found_error = True
                            print('#### error\n%s\n' %  run_data['error'])
                        if not found_error:
                            print('No errors found.')
                            print("Run 'uap %s fix-problems --first-error' to investigate.'"
                                    % p.args.config.name)
                            print('')
            else:
                print("Some tasks are bad. Run 'uap %s status --details' to see the details." %
                        p.args.config.name)
        # now check ping files and print some warnings and instructions if
        # something's fishy
        p.check_ping_files(print_more_warnings = True if args.verbose > 0 else False)

        # Now check whether we can volatilize files, but don't do it.
        p.check_volatile_files()

