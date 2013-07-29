#!./python_env/bin/python

import sys
sys.path.append('./include')
import abstract_step
import argparse
import os
import pipeline
import subprocess
import yaml

'''
This script uses graphviz to produce graphs that display information about the 
tasks processed by the pipeline. 
'''

parser = argparse.ArgumentParser(
    description="This script displays by default information about all tasks " +
                "of the pipeline as configured in 'config.yaml'. But the " +
                "displayed information can be narrowed down via command " +
                "line options.",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("--all",
                    dest="all",
                    action="store_true",
                    default=False,
                    help="Renders all ")

parser.add_argument("-t","--task",
                    dest="task",
                    nargs='*',
                    default=list(),
                    type=str,
                    help="Displays only the named task IDs" +
                    "Can take multiple task ID(s) as input. A task ID " +
                    "looks like ths 'step_name/run_id'. A list of all " +
                    "task IDs is returned by running './status.py'.")

args = parser.parse_args()

def escape(s):
    result = ''
    for c in s:
        result += "x%x" % ord(c)
    return result

GRADIENTS = {
    'burn': [
        [0.0, '#ffffff'],
        [0.2, '#fce94f'],
        [0.4, '#fcaf3e'],
        [0.7, '#a40000'],
        [1.0, '#000000']
    ],
    'green': [
        [0.0, '#ffffff'],
        [1.0, '#4e9a06']
    ],
    'traffic_lights': [
        [0.0, '#d5291a'],
        [0.5, '#fce94f'],
        [1.0, '#73a946']
    ]
}


def mix(a, b, amount):
    rA = float(int(a[1:3], 16)) / 255.0
    gA = float(int(a[3:5], 16)) / 255.0
    bA = float(int(a[5:7], 16)) / 255.0
    rB = float(int(b[1:3], 16)) / 255.0
    gB = float(int(b[3:5], 16)) / 255.0
    bB = float(int(b[5:7], 16)) / 255.0
    rC = rB * amount + rA * (1.0 - amount)
    gC = gB * amount + gA * (1.0 - amount)
    bC = bB * amount + bA * (1.0 - amount)
    result = '#%02x%02x%02x' % (int(rC * 255.0), int(gC * 255.0), int(bC * 255.0))
    return result


def gradient(x, gradient):
    x = max(x, 0.0)
    x = min(x, 1.0)
    i = 0
    while (i < len(gradient) - 2 and gradient[i + 1][0] < x):
        i += 1
    colorA = gradient[i][1]
    colorB = gradient[i + 1][1]
    return mix(colorA, colorB, (x - gradient[i][0]) / (gradient[i + 1][0] - gradient[i][0]))


def main():
    p = pipeline.Pipeline(arguments=args)
    
    if args.all or args.task:
        logs = []
        if args.all:
            for task in p.task_for_task_id.values():
                annotation_path = os.path.join(task.step.get_output_directory(), '.%s-annotation.yaml' % task.run_id)
                if os.path.exists(annotation_path):
                    log = yaml.load(open(annotation_path))
                    logs.append(log)
        else:
            for task_id in args.task:
                task = p.task_for_task_id[task_id]
                annotation_path = os.path.join(task.step.get_output_directory(), '.%s-annotation.yaml' % task.run_id)
                if os.path.exists(annotation_path):
                    log = yaml.load(open(annotation_path))
                    logs.append(log)
                else:
                    print("Unable to find annotation at %s." % annotation_path)
        gv = abstract_step.AbstractStep.render_pipeline(logs)
        with open('out.gv', 'w') as f:
            f.write(gv)
            
        exit(0)
    
    dot = subprocess.Popen(['dot', '-Tsvg'], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    
    f = dot.stdin
    
    f.write("digraph {\n")
    f.write("  rankdir = TB;\n")
    f.write("  splines = true;\n")
    f.write("    graph [fontname = Helvetica, fontsize = 12, size = \"14, 11\", nodesep = 0.2, ranksep = 0.3];\n")
    f.write("    node [fontname = Helvetica, fontsize = 12, shape = rect];\n")
    f.write("    edge [fontname = Helvetica, fontsize = 12];\n")
    for step_name, step in p.steps.items():
        total_runs = len(step.get_run_ids())
        finished_runs = 0
        for _ in step.get_run_ids():
            if step.get_run_state(_) == p.states.FINISHED:
                finished_runs += 1

        f.write("subgraph cluster_%s {\n" % step_name)
        
        label = step_name
        if step_name != step.__module__:
            label = "%s\\n(%s)" % (step_name, step.__module__)
        f.write("    %s [label=\"%s\", style = filled, fillcolor = \"#fce94f\"];\n" % (step_name, label))
        color = gradient(float(finished_runs) / total_runs if total_runs > 0 else 0.0, GRADIENTS['traffic_lights'])
        color = mix(color, '#ffffff', 0.5)
        f.write("    %s_progress [label=\"%1.0f%%\", style = filled, fillcolor = \"%s\" height = 0.3];\n" % (step_name, float(finished_runs) * 100.0 / total_runs if total_runs > 0 else 0.0, color))
        f.write("    %s -> %s_progress [arrowsize = 0];\n" % (step_name, step_name))
        f.write("    {rank=same; %s %s_progress}\n" % (step_name, step_name))
        
        for c in step._connections:
            connection_key = escape(('%s/%s' % (step_name, c)).replace('/', '__'))
            f.write("    %s [label=\"%s\", shape = ellipse, fontsize = 10];\n" % (connection_key, c))
            if c[0:3] == 'in/':
                f.write("    %s -> %s;\n" % (connection_key, step_name))
            else:
                f.write("    %s -> %s;\n" % (step_name, connection_key))
                
        f.write("  graph[style=dashed];\n")
        f.write("}\n")
            
    for step_name, step in p.steps.items():
        for other_step in step.dependencies:
            #f.write("    %s -> %s;\n" % (other_step.get_step_name(), step_name))
            
            for in_key in step._connections:
                if in_key[0:3] != 'in/':
                    continue
                
                out_key = in_key.replace('in/', 'out/')
                allowed_steps = None
                if '_connect' in step.options:
                    if in_key in step.options['_connect']:
                        declaration = step.options['_connect'][in_key]
                        if declaration.__class__ == str:
                            if '/' in declaration:
                                parts = declaration.split('/')
                                allowed_steps = set()
                                allowed_steps.add(parts[0])
                                out_key = 'out/' + parts[1]
                            else:
                                out_key = 'out/' + declaration
                        else:
                            raise StandardError("Invalid _connect value: %s" % yaml.dump(declaration))
                        
                for real_outkey in other_step._connections:
                    if real_outkey[0:4] != 'out/':
                        continue
                    if out_key == real_outkey:
                        connection_key = escape(('%s/%s' % (step_name, in_key)).replace('/', '__'))
                        other_connection_key = escape(('%s/%s' % (other_step.get_step_name(), out_key)).replace('/', '__'))
                        f.write("    %s -> %s;\n" % (other_connection_key, connection_key))
    f.write("}\n")
    
    dot.stdin.close()
    
    with open('steps.svg', 'w') as f:
        f.write(dot.stdout.read())
        
if __name__ == '__main__':
    main()
