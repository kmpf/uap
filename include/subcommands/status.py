#!/usr/bin/env python
# encoding: utf-8

import sys
from contextlib import closing
import logging
import pydoc
import string
from cStringIO  import StringIO
import yaml
from deepdiff import DeepDiff
from tqdm import tqdm

import pipeline
from uaperrors import UAPError

'''
By default, this script displays information about all tasks of the pipeline
configured in 'config.yaml'. But the displayed information can be narrowed 
down via command line options.

'''

logger = logging.getLogger("uap_logger")

def main(args):
    p = pipeline.Pipeline(arguments=args)
    group_by_status = True

    if args.sources:
        # print all sources (i. e. instances of AbstractSourceStep)
        p.print_source_runs()

    elif args.job_ids:
        ids = p.get_cluster_job_ids()
        print(' '.join(ids))
        return

    elif len( args.run ) >= 1:
        # print run infos of one or more specific tasks
        for task_id in args.run:
            parts = task_id.split('/')
            if len(parts) != 2:
                raise StandardError("Invalid task ID %s." % task_id)
            step_name = parts[0]
            run_id = parts[1]
            report = p.steps[step_name].get_run(run_id).as_dict()
            report['state'] = p.steps[step_name].get_run_state(run_id)
            print(yaml.dump(report, default_flow_style = False))
        
    elif args.graph:
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

            line = "%s%s%s [%s]" % (''.join(_).replace("─└", "─┴"), step.get_step_name(), original_step_name_label, step.get_run_info_str())
            print(line)
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
        output = list()
        tasks_for_status = {}
        tasks = p.all_tasks_topologically_sorted
        for task in tqdm(tasks, desc='tasks'):
            state = task.get_task_state()
            if not state in tasks_for_status:
                tasks_for_status[state] = list()
            tasks_for_status[state].append(task)
            if not group_by_status:
                output.append(
                    "[%s] %s" % (state[0].lower(), task))
                #print("[%s] %s" % (task.get_task_state()[0].lower(), task))
        if group_by_status:
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

        if p.states.CHANGED in tasks_for_status.keys():
            if args.details:
                print('')
                for task in tasks_for_status[p.states.CHANGED]:
                    heading = 'changes in task %s' % task
                    print(heading)
                    print('-'*len(heading))
                    run = task.get_run()
                    anno_file = run.get_annotation_path()
                    try:
                        with open(anno_file, 'r') as fl:
                            anno_data = yaml.load(fl, Loader=yaml.FullLoader)
                    except IOError:
                        print('The annotation file could not be read: %s.' %
                                anno_file)
                    else:
                        old_strcut = anno_data['run']['structure']
                        new_struct = run.get_run_structure()
                        diff =  DeepDiff(old_strcut, new_struct)
                        print(yaml.dump(dict(diff)))
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
                    try:
                        with open(anno_file, 'r') as fl:
                            anno_data = yaml.load(fl, Loader=yaml.FullLoader)
                    except IOError:
                        print('The annotation file could not be read: %s.\n' %
                                anno_file)
                    else:
                        failed = dict()
                        try:
                            procs = anno_data['pipeline_log']['processes']
                            for proc in procs:
                                if proc.get('exit_code', 0) == 0:
                                    continue
                                failed[proc['name']] = proc['stderr_copy']['tail']
                        except KeyError as e:
                            print('The annotation file "%s" seems badly '
                                    'formated: %s\n' % (anno_file, e))
                        if failed:
                            found_error = True
                            print('stderr of failed commands:')
                            print(yaml.dump(failed))
                        else:
                            print('No failed commands found in the annotation file.\n')
                        if 'error' in anno_data:
                            found_error = True
                            print('uap error:\n%s\n' % anno_data['error'])
                        if not found_error:
                            print('No errors found.')
                            print("Run 'uap %s fix-problems --first-error' to investigate.'"
                                    % p.args.config.name)
                            print('')
            else:
                print("Some tasks are bad. Run 'uap %s status --details' to see the details." %
                        p.args.config.name)
            print("If you want to force overwrite of the bad tasks, run\n"
                  "'uap %s run-locally --force' or 'uap %s submit-to-cluster --force'." %
                  (p.args.config.name, p.args.config.name))
    # now check ping files and print some warnings and instructions if
    # something's fishy
    p.check_ping_files(print_more_warnings = True if args.verbose > 0 else False)
    
    # Now check whether we can volatilize files, but don't do it.
    p.check_volatile_files()
    
