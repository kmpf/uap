#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import subprocess
import yaml

def main():
    p = pipeline.Pipeline()
    
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

        label = step_name
        if step_name != step.__module__:
            label = "%s\\n(%s)" % (step_name, step.__module__)
        f.write("    %s [label=\"%s\", style = filled, fillcolor = \"#fce94f\"];\n" % (step_name, label))
        f.write("    %s_progress [label=\"%1.1f%%\", style = filled, height = 0.3];\n" % (step_name, float(finished_runs) * 100.0 / total_runs))
        f.write("    %s -> %s_progress [arrowsize = 0];\n" % (step_name, step_name))
        f.write("    {rank=same; %s %s_progress}\n" % (step_name, step_name))
        
        for c in step._connections:
            connection_key = ('%s/%s' % (step_name, c)).replace('/', '__')
            f.write("    %s [label=\"%s\", shape = ellipse, fontsize = 10];\n" % (connection_key, c))
            if c[0:3] == 'in/':
                f.write("    %s -> %s [style=dashed, color = \"#888a85\"];\n" % (connection_key, step_name))
            else:
                f.write("    %s -> %s [style=dashed, color = \"#888a85\"];\n" % (step_name, connection_key))
            
    for step_name, step in p.steps.items():
        for other_step in step.dependencies:
            f.write("    %s -> %s;\n" % (other_step.get_step_name(), step_name))
            for c in step._connections:
                if c[0:3] != 'in/':
                    continue
                connection_key = ('%s/%s' % (step_name, c)).replace('/', '__')
                for oc in other_step._connections:
                    if oc[0:4] != 'out/':
                        continue
                    other_connection_key = ('%s/%s' % (other_step.get_step_name(), oc)).replace('/', '__')
                    if c.replace('in/', 'out/') == oc:
                        f.write("    %s -> %s [style=dashed, color = \"#888a85\"];\n" % (other_connection_key, connection_key))
    f.write("}\n")
    
    dot.stdin.close()
    
    with open('steps.svg', 'w') as f:
        f.write(dot.stdout.read())
        
if __name__ == '__main__':
    main()
