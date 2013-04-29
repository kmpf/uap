#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import yaml

def main():
    p = pipeline.Pipeline()

    if len(sys.argv) > 1:
        if sys.argv[1] == '--sources':
            # print all sources (i. e. instances of AbstractSourceStep)
            p.print_source_runs()
        else:
            # print a specific task
            for task_id in sys.argv[1:]:
                parts = task_id.split('/')
                if len(parts) != 2:
                    raise StandardError("Invalid task ID %s." % task_id)
                step_name = parts[0]
                run_id = parts[1]
                print(yaml.dump(p.steps[step_name].get_run_info()[run_id], default_flow_style = False))
    else:
        # print all tasks
        p.print_tasks()
        
if __name__ == '__main__':
    main()
