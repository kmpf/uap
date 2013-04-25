#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import yaml

def main():
    p = pipeline.Pipeline()

    if len(sys.argv) > 1:
        if sys.argv[1] == '--samples':
            # print all samples
            print yaml.dump(p.all_samples, default_flow_style = False)
            print yaml.dump(p.steps[0].get_run_info(), default_flow_style = False)
        else:
            # print a specific task
            task = p.task_for_task_id[sys.argv[1]]
            print(yaml.dump(task.step.get_run_info()[task.run_id], default_flow_style = False))
    else:
        # print all tasks
        p.print_tasks()
        
if __name__ == '__main__':
    main()
