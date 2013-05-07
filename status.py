#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import string
import yaml

def main():
    p = pipeline.Pipeline()
    
    group_by_status = False
    if '--group' in sys.argv:
        group_by_status = True
        sys.argv.remove('--group')
    
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
        # print all tasks
        '''
        prints a summary of all tasks, indicating whether each taks is
          - ``[r]eady``
          - ``[w]aiting``
          - ``[f]inished``
        '''
        tasks_for_status = {}
        for task in p.all_tasks:
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
                for task in tasks_for_status[status]:
                    print("[%s] %s" % (task.get_task_state()[0].lower(), task))
                print('')
        print("tasks: %d total, %s" % (len(p.all_tasks), ', '.join(["%d %s" % (len(tasks_for_status[_]), _.lower()) for _ in p.states.order if _ in tasks_for_status])))
        
if __name__ == '__main__':
    main()