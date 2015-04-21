#!/usr/bin/env python
# encoding: utf-8

import sys
from contextlib import closing
import logging
import pydoc
import string
from cStringIO  import StringIO
import yaml

import pipeline

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

    elif len( args.task ) >= 1:
        # print run infos of one or more specific tasks
        for task_id in args.task:
            parts = task_id.split('/')
            if len(parts) != 2:
                raise StandardError("Invalid task ID %s." % task_id)
            step_name = parts[0]
            run_id = parts[1]
            report = p.steps[step_name].get_run_info()[run_id].as_dict()
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

            line = "%s%s%s [%s]" % (''.join(_).replace("─└", "─┴"), step.get_step_name(), original_step_name_label, step.get_run_info_str())
            print(line)
    else:
        # print all tasks
        '''
        prints a summary of all tasks, indicating whether each taks is
        - ``[r]eady``
        - ``[w]aiting``
        - ``[q]ueued``
        - ``[e]xecuting``
        - ``[f]inished``
        '''
        output = list()
        tasks_for_status = {}
        for task in p.all_tasks_topologically_sorted:
            state = task.get_task_state()
            if not state in tasks_for_status:
                tasks_for_status[state] = list()
            tasks_for_status[state].append(task)
            if not group_by_status:
                output.append(
                    "[%s] %s" % (task.get_task_state()[0].lower(), task))
                #print("[%s] %s" % (task.get_task_state()[0].lower(), task))
        if group_by_status:
            for status in p.states.order:
                if not status in tasks_for_status:
                    continue
                heading = "%s tasks" % string.capwords(status)
                output.append(heading)
                output.append('-' * len(heading))
                if args.summarize:
                    step_count = dict()
                    step_order = list()
                    for task in tasks_for_status[status]:
                        if not str(task.step) in step_count:
                            step_count[str(task.step)] = 0
                            step_order.write(str(task.step))
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
                                         task.get_task_state()[0].lower(),
                                         task))
                    output.append('')
            output.append("tasks: %d total, %s"
                          % (len(p.all_tasks_topologically_sorted),
                             ', '.join(["%d %s" % (
                                 len(tasks_for_status[_]),
                                 _.lower()) for _ in p.states.order \
                                        if _ in tasks_for_status])))
            pydoc.pager("\n".join(output))
    # now check ping files and print some warnings and instructions if
    # something's fishy
    p.check_ping_files(print_more_warnings = True if args.verbose > 0 else False)
    
    # Now check whether we can volatilize files, but don't do it.
    p.check_volatile_files()
    
