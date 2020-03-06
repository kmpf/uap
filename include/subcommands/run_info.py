#!/usr/bin/env python
# encoding: utf-8

import sys
import logging
import yaml
from subprocess import list2cmdline

import pipeline
import pipeline_info
import command as command_info

'''
By default, this script displays information about all runs of the pipeline
configured in 'config.yaml'. But the displayed information can be narrowed 
down via command line options.

'''

logger = logging.getLogger("uap_logger")

def main(args):
    args.no_tool_checks = True
    p = pipeline.Pipeline(arguments=args)
    group_by_status = True

    if args.sources:
        # print all sources (i. e. instances of AbstractSourceStep)
        p.print_source_runs()

    else:
        # print run infos of one or more specific tasks
        for task_id in p.get_task_with_list(as_string=True, exclusive=True):
            parts = task_id.split('/')
            if len(parts) != 2:
                raise StandardError("Invalid run ID %s." % task_id)
            step_name = parts[0]
            run_id = parts[1]
            run = p.steps[step_name].get_run(run_id)
            report = run.as_dict()
            shebang = "#!/usr/bin/env bash"
            print(shebang + "\n")
            report_header = "%s/%s -- Report" % (step_name, run_id)
            print("# " + report_header)
            print("# " + "=" * len(report_header) + "\n#")
            dump = yaml.dump(report, default_flow_style = False)
            for line in dump.split('\n'):
                print("# " + line)
            exec_header = "%s/%s -- Commands" % (step_name, run_id)
            print("# " + exec_header)
            print("# " + "=" * len(exec_header) + "\n")
            eg_count = 1
            for exec_group in run.get_exec_groups():
                goc_header = "%d. Group of Commands" % eg_count
                eg_count += 1
                pipe_count = 1
                cmd_count = 1
                for poc in exec_group.get_pipes_and_commands():
                    # for each pipe or command (poc)
                    # check if it is a pipeline ...
                    
                    if isinstance(poc, pipeline_info.PipelineInfo):
                        cmd_header = goc_header + " -- %d. Pipeline" % pipe_count
                        pipe_count += 1
                    elif isinstance(poc, command_info.CommandInfo):
                        cmd_header = goc_header + " -- %d. Command" % cmd_count
                        cmd_count += 1
                    print("# " + cmd_header)
                    print("# " + "-" * len(cmd_header) + "\n")
                    print(poc.get_command_string() + "\n")
