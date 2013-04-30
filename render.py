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

        f.write("subgraph cluster_%s {\n" % step_name)
        
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
                
        f.write("  graph[style=dashed];\n")
        f.write("}\n")
            
    for step_name, step in p.steps.items():
        for other_step in step.dependencies:
            f.write("    %s -> %s;\n" % (other_step.get_step_name(), step_name))
            
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
                        connection_key = ('%s/%s' % (step_name, in_key)).replace('/', '__')
                        other_connection_key = ('%s/%s' % (other_step.get_step_name(), out_key)).replace('/', '__')
                        f.write("    %s -> %s [style=dashed, color = \"#888a85\"];\n" % (other_connection_key, connection_key))
    f.write("}\n")
    
    dot.stdin.close()
    
    with open('steps.svg', 'w') as f:
        f.write(dot.stdout.read())
        
if __name__ == '__main__':
    main()
