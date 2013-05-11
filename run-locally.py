#!./python_env/bin/python

import sys
sys.path.append('./include')
import copy
import pipeline
import process_pool
import yaml


def main():
    p = pipeline.Pipeline()

    if len(sys.argv) > 1:
        # execute the specified tasks
        for task_id in sys.argv[1:]:
            task = p.task_for_task_id[task_id]
            if task.get_task_state() == p.states.FINISHED:
                continue
            if task.get_task_state() in [p.states.READY, p.states.QUEUED]:
                task.run()
            else:
                raise StandardError("Unexpected task state for %s: %s" % (task, task.get_task_state()))
    else:
        # execute all tasks
        for task in p.all_tasks_topologically_sorted:
            if task.get_task_state() == p.states.FINISHED:
                continue
            if task.get_task_state() in [p.states.READY, p.states.QUEUED]:
                task.run()
            else:
                raise StandardError("Unexpected task state for %s: %s" % (task, task.get_task_state()))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Script terminated by user.")
    finally:
        # make sure all child processes get terminated
        process_pool.ProcessPool.kill_all_child_processes()
