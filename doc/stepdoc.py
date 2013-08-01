#!../python_env/bin/python

import sys
sys.path.append('../include')
sys.path.append('../include/steps')
sys.path.append('../include/sources')
import abstract_step
import glob
import os
import pipeline
import string
import yaml

def doc_module(module_name, fout):
    step_class = abstract_step.AbstractStep.get_step_class_for_key(module_name)
    step = step_class(None)
    fout.write(module_name + "\n")
    fout.write('~' * len(module_name) + "\n\n")
    if step.__doc__:
        doc = step.__doc__.split("\n")
        for line in doc:
            fout.write(line.strip() + "\n")
        
    # print connections
    fout.write("**Connections:**\n\n")
    fout.write(".. graphviz::\n")
    fout.write("\n")
    fout.write("   digraph foo {\n")
    fout.write("      rankdir = LR;\n")
    fout.write("      splines = true;\n")
    fout.write("      graph [fontname = Helvetica, fontsize = 12, size = \"14, 11\", nodesep = 0.2, ranksep = 0.3];\n")
    fout.write("      node [fontname = Helvetica, fontsize = 12, shape = rect];\n")
    fout.write("      edge [fontname = Helvetica, fontsize = 12];\n")
    fout.write("      %s [style=filled, fillcolor=\"#fce94f\"];\n" % module_name)
    for index, c in enumerate(sorted(step._connections)):
        c = c.split('/')
        if c[0] == 'out':
            fout.write("      out_%d [label=\"%s\"];\n" % (index, c[1]))
            fout.write("      %s -> out_%d;\n" % (module_name, index))
        else:
            fout.write("      in_%d [label=\"%s\"];\n" % (index, c[1]))
            fout.write("      in_%d -> %s;\n" % (index, module_name))
    fout.write("   }    \n")
    fout.write("\n")
    
    # print options
    fout.write("**Options:**\n")
    for key in sorted(step._defined_options.keys()):
        option = step._defined_options[key]
        print(option)
        fout.write("  - **%s** (%s, %s)" % (
            key, 
            '/'.join([_.__name__ for _ in option['types']]),
            'optional' if option['optional'] else 'required'
            ))
        if option['description']:
            fout.write(" -- %s" % option['description'])
        fout.write("\n")
        fout.write("    \n")
        if option['choices']:
            fout.write("    - possible values:\n")
            fout.write("    \n")
            for v in sorted(option['choices']):
                fout.write("      - %s\n" % v)
            
    fout.write("\n")

    # print tools
    if len(step._tools) > 0:
        fout.write("**Required tools:** %s\n" % ', '.join(sorted(step._tools.keys())))
        fout.write("\n")
        
    if not abstract_step.AbstractSourceStep in step.__class__.__bases__:
        # this is not a source step, print cores
        fout.write("**Cores:** %s\n" % step._cores)
        fout.write("\n")
        
    
    '''
    print("Cores: %d" % step._cores)
    print("Connections: %s" % step._connections)
    print("Connection restrictions: %s" % step._connection_restrictions)
    print("Tools: %s" % step._tools)
    print("Options: %s" % sorted(step._defined_options.keys()))
    print(step.__doc__)
    '''

def main():
    with open('source/steps.rst', 'w') as fout:
        fout.write("Available steps\n")
        fout.write("===============\n\n")

        fout.write("Source steps\n")
        fout.write("------------\n\n")
        modules = glob.glob('../include/sources/*.py')
        for m in sorted(modules):
            module_name = os.path.basename(m).replace('.py', '')
            doc_module(module_name, fout)

        fout.write("Processing steps\n")
        fout.write("----------------\n\n")
        modules = glob.glob('../include/steps/*.py')
        for m in sorted(modules):
            module_name = os.path.basename(m).replace('.py', '')
            if module_name == 'io_step':
                continue
            doc_module(module_name, fout)

if __name__ == '__main__':
    main()
