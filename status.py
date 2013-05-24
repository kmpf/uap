#!./python_env/bin/python
# encoding: utf-8

import sys
sys.path.append('./include')
import pipeline
import string
import yaml

def main():
    p = pipeline.Pipeline()
    
    group_by_status = True
    summarize = False
    graph = False
    
    if '--summarize' in sys.argv:
        sys.argv.remove('--summarize')
        summarize = True

    if '--graph' in sys.argv:
        sys.argv.remove('--graph')
        graph = True
    
    if len(sys.argv) > 1:
        if sys.argv[1] == '--sources':
            # print all sources (i. e. instances of AbstractSourceStep)
            p.print_source_runs()
        else:
            # print one or more specific tasks
            for task_id in sys.argv[1:]:
                parts = task_id.split('/')
                if len(parts) != 2:
                    raise StandardError("Invalid task ID %s." % task_id)
                step_name = parts[0]
                run_id = parts[1]
                report = p.steps[step_name].get_run_info()[run_id]
                report['state'] = p.steps[step_name].get_run_state(run_id)
                print(yaml.dump(report, default_flow_style = False))
    else:
        if graph:
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

            lines = ["%s%s (%s)" % (''.join(_).replace("─└", "─┴"), step_order[index], p.steps[step_order[index]].get_run_info_str()) for index, _ in enumerate(lines)]
            for line in lines:
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
            tasks_for_status = {}
            for task in p.all_tasks_topologically_sorted:
                state = task.get_task_state()
                if not state in tasks_for_status:
                    tasks_for_status[state] = list()
                tasks_for_status[state].append(task)
                if not group_by_status:
                    print("[%s] %s" % (task.get_task_state()[0].lower(), task))
            if group_by_status:
                for status in p.states.order:
                    if not status in tasks_for_status:
                        continue
                    heading = "%s tasks" % string.capwords(status)
                    print(heading)
                    print('-' * len(heading))
                    if summarize:
                        step_count = dict()
                        step_order = list()
                        for task in tasks_for_status[status]:
                            if not str(task.step) in step_count:
                                step_count[str(task.step)] = 0
                                step_order.append(str(task.step))
                            step_count[str(task.step)] += 1
                        for step_name in step_order:
                            print("[%s]%4d %s" % (status.lower()[0], step_count[step_name], step_name))
                    else:
                        for task in tasks_for_status[status]:
                            print("[%s] %s" % (task.get_task_state()[0].lower(), task))
                    print('')
            print("tasks: %d total, %s" % (len(p.all_tasks_topologically_sorted), ', '.join(["%d %s" % (len(tasks_for_status[_]), _.lower()) for _ in p.states.order if _ in tasks_for_status])))
            
    # now check ping files and print some warnings and instructions if something's fishy
    p.check_ping_files()
    
    # Now check whether we can volatilize files.
    # This happens when you declare a step as volatile after all downstream 
    # steps have been run - without the _volatile declaration, the files will
    # remain untouched. 
    p.check_volatile_files()
        
if __name__ == '__main__':
    main()