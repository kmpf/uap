#!/usr/bin/env python
# encoding: utf-8

import sys
import copy
import glob
import logging
import os
import re
import socket
import StringIO
import subprocess
import textwrap
import yaml

import abstract_step
import pipeline
import misc
'''
This script uses graphviz to produce graphs that display information about the 
tasks processed by the pipeline. 
'''

logger = logging.getLogger("uap_logger")

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


def main(args):
    p = pipeline.Pipeline(arguments=args)

    # Test if dot is available
    dot_version = ['dot', '-V']
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(dot_version, stdout = devnull)
    except subprocess.CalledProcessError as e:
        raise StandardError("Execution of %s failed. GraphViz seems to be "
                            "unavailable." % " ".join(dot_version))

    if args.files:
        logger.info("Going to plot the graph containing all files of the analysis")
        raise StandardError("Sorry, feature not implemented yet!")
    elif args.steps:
        logger.info("Create a graph showing the DAG of the analysis")

        render_graph_for_all_steps(p, args)

    elif args.step_task:
        if len(args.step_task) >= 1:
            # Compile a list of all tasks
            task_list = list()
            for task_id in args.step_task:
                if '/' in task_id:
                    task = p.task_for_task_id[task_id]
                    task_list.append(task)
                else:
                    for task in p.all_tasks_topologically_sorted:
                        if str(task)[0:len(task_id)] == task_id:
                            task_list.append(task)

            for task in task_list:
                outdir = task.get_run().get_output_directory()
                anno_files = glob.glob(os.path.join(
                    outdir, ".%s*.annotation.yaml" % task.get_run().get_run_id()
                ))

                for f in anno_files:
                    logger.info("Going to plot the graph for task: %s" % task)
                    render_single_annotation(f)

    else:
        # Just find everything that looks like an annotation file and generate
        # the plots for it.
        directory = os.path.abspath(p.config['destination_path'])
        if not os.path.isdir(directory):
            raise StandardError("%s is not a directory." % directory)

        annotation_files = list()
        for root, dirs, files in os.walk(directory):
            anno_file_regex = re.compile('^\..*-annotation-.{6}\.yaml$')
            for f in files:
                if re.search(anno_file_regex, f):
                    annotation_files.append(os.path.join(root, f))

        for annotation_file in annotation_files:
            render_single_annotation(annotation_file)
        


def render_graph_for_all_steps(p, args):
    configuration_path = p.get_config_filepath()
    if args.simple:
        svg_file = configuration_path.replace('.yaml', '.simple.svg')
        png_file = configuration_path.replace('.yaml', '.simple.png')        
    else:
        svg_file = configuration_path.replace('.yaml', '.svg')
        png_file = configuration_path.replace('.yaml', '.png')

    dot = subprocess.Popen(['dot', '-Tsvg'],
                           stdin = subprocess.PIPE,
                           stdout = subprocess.PIPE)
    
    f = dot.stdin
    
    f.write("digraph {\n")
    f.write("  rankdir = TB;\n")
    f.write("  splines = true;\n")
    f.write("    graph [fontname = Helvetica, fontsize = 12, size = \"14, 11\", "
            "nodesep = 0.2, ranksep = 0.3];\n")
    f.write("    node [fontname = Helvetica, fontsize = 12, shape = rect];\n")
    f.write("    edge [fontname = Helvetica, fontsize = 12];\n")
    for step_name, step in p.get_steps().items():
        total_runs = len(step.get_run_ids())
        finished_runs = 0
        for _ in step.get_run_ids():
            if step.get_run_state(_) == p.states.FINISHED:
                finished_runs += 1

        f.write("subgraph cluster_%s {\n" % step_name)
        
        label = step_name
        if step_name != step.__module__:
            label = "%s\\n(%s)" % (step_name, step.__module__)
        f.write("    %s [label=\"%s\", style = filled, fillcolor = \"#fce94f\"];\n"
                % (step_name, label))
        color = gradient(float(finished_runs) / total_runs \
                         if total_runs > 0 else 0.0, GRADIENTS['traffic_lights'])
        color = mix(color, '#ffffff', 0.5)
        f.write("    %s_progress [label=\"%s/%s\", style = filled, "
                "fillcolor = \"%s\" height = 0.3];\n"
                % (step_name, finished_runs, total_runs, color))
        f.write("    %s -> %s_progress [arrowsize = 0];\n"
                % (step_name, step_name))
        f.write("    {rank=same; %s %s_progress}\n" % (step_name, step_name))

        if not args.simple:
            for c in step._connections:
                connection_key = escape(('%s/%s'
                                         % (step_name, c)).replace('/', '__'))
                f.write("    %s [label=\"%s\", shape = ellipse, fontsize = 10];\n"
                        % (connection_key, c))
                if c[0:3] == 'in/':
                    f.write("    %s -> %s;\n" % (connection_key, step_name))
                else:
                    f.write("    %s -> %s;\n" % (step_name, connection_key))
                
        f.write("  graph[style=dashed];\n")
        f.write("}\n")
            
    for step_name, step in p.steps.items():
        for other_step in step.dependencies:
            if args.simple:
                f.write("    %s -> %s;\n"
                        % (other_step.get_step_name(), step_name))
            else:
            
                for in_key in step._connections:
                    if in_key[0:3] != 'in/':
                        continue
                
                    out_key = in_key.replace('in/', 'out/')
                    allowed_steps = None
                    if '_connect' in step.get_options():
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
                                raise StandardError("Invalid _connect value: %s"
                                                    % yaml.dump(declaration))
                        
                    for real_outkey in other_step._connections:
                        if real_outkey[0:4] != 'out/':
                            continue
                        if out_key == real_outkey:
                            connection_key = escape(
                                ('%s/%s' % (step_name, in_key)).replace('/', '__')
                            )
                            other_connection_key = escape(
                                ('%s/%s' % (other_step.get_step_name(),
                                            out_key)).replace('/', '__')
                            )
                            f.write("    %s -> %s;\n"
                                    % (other_connection_key, connection_key))

    f.write("}\n")
    
    dot.stdin.close()

    svg = dot.stdout.read()
    with open(svg_file, 'w') as f:
        f.write(svg)


    #print(yaml.dump(pipeline.get_steps()))

    

def render_single_annotation(annotation_path):
    dot_file = annotation_path.replace('.yaml', '.dot')
    # Replace leading dot to make graphs easier to find
    annotation_path = annotation_path[0].replace('.', '')
    svg_file = annotation_path.replace('.yaml', '.svg')
    png_file = annotation_path.replace('.yaml', '.png')

    log = dict()
    with open(annotation_path, 'r') as f:
        log = yaml.load(f)
    try:
        gv = render_pipeline([log])
        with open(dot_file, 'w') as f:
            f.write(gv)
        
        dot = subprocess.Popen(['dot', '-Tsvg', '-o%s' % svg_file, dot_file],
                               stdin = subprocess.PIPE,
                               stdout = subprocess.PIPE)
#        dot.stdin.write(gv)
#        dot.stdin.close()
#        svg = dot.stdout.read()
#        with open(svg_file, 'w') as f:
#            f.write(svg)
            
        dot = subprocess.Popen(['dot', '-Tpng', '-o%s' % png_file, dot_file],
                               stdin = subprocess.PIPE,
                               stdout = subprocess.PIPE)
#        dot.stdin.write(gv)
#        dot.stdin.close()
#        png = dot.stdout.read()
#        with open(png_file, 'w') as f:
#            f.write(png)
    except:
        print(sys.exc_info())
        import traceback
        traceback.print_tb(sys.exc_info()[2])
        pass


def render_pipeline(logs):
    hash = {'nodes': {}, 'edges': {}, 'clusters': {}, 'graph_labels': {}}
    for log in logs:
        temp = render_pipeline_hash(log)
        for _ in ['nodes', 'edges', 'clusters', 'graph_labels']:
            hash[_].update(temp[_])

    f = StringIO.StringIO()
    f.write("digraph {\n")
    f.write("    rankdir = TB;\n")
    f.write("    splines = true;\n")
    f.write("    graph [fontname = Helvetica, fontsize = 12, size = "
            "\"14, 11\", nodesep = 0.2, ranksep = 0.3, labelloc = t, "
            "labeljust = l];\n")
    f.write("    node [fontname = Helvetica, fontsize = 12, shape = rect, "
            "style = filled];\n")
    f.write("    edge [fontname = Helvetica, fontsize = 12];\n")
    f.write("\n")
    
    f.write("    // nodes\n")
    f.write("\n")
    for node_key, node_info in hash['nodes'].items():
        f.write("    _%s" % node_key)
        if len(node_info) > 0:
            f.write(" [%s]" % ', '.join(['%s = "%s"' % (k, node_info[k]) \
                                         for k in node_info.keys()]))
        f.write(";\n")
            
    f.write("\n")
        
    f.write("    // edges\n")
    f.write("\n")
    for edge_pair in hash['edges'].keys():
        if edge_pair[0] in hash['nodes'] and edge_pair[1] in hash['nodes']:
            f.write("    _%s -> _%s;\n" % (edge_pair[0], edge_pair[1]))
        
    f.write("\n")
        
        
    if len(hash['graph_labels']) == 1:
        f.write("    graph [label=\"%s\"];\n" % 
                hash['graph_labels'].values()[0])
    f.write("}\n")
        
    result = f.getvalue()
    f.close()
    return result
        
def render_pipeline_hash(log):
        
    def pid_hash(pid, suffix = ''):
        hashtag = "%s/%s/%d/%s" % (log['step']['name'], 
                                   log['run']['run_id'], 
                                   pid, suffix)
        return misc.str_to_sha1(hashtag)
        
    def file_hash(path):
        if path in log['step']['known_paths']:
            if 'real_path' in log['step']['known_paths'][path]:
                path = log['step']['known_paths'][path]['real_path']
        return misc.str_to_sha1(path)
        
        
    pipe_hash = dict()
    pipe_hash['nodes'] = dict()
    pipe_hash['edges'] = dict()
    pipe_hash['clusters'] = dict()
    pipe_hash['graph_labels'] = dict()
        
    def add_file_node(path):
        if not path in log['step']['known_paths']:
            return
                
        if 'real_path' in log['step']['known_paths'][path]:
            path = log['step']['known_paths'][path]['real_path']
        #label = log['step']['known_paths'][path]['label']
        color = '#ffffff'
        if log['step']['known_paths'][path]['type'] == 'fifo':
            color = '#c4f099'
            label = os.path.basename(path)
        elif log['step']['known_paths'][path]['type'] == 'file':
            color = '#8ae234'
            label = os.path.basename(path)
        elif log['step']['known_paths'][path]['type'] == 'step_file':
            color = '#97b7c8'
            label = log['step']['known_paths'][path]['label']
            if path in log['step']['known_paths']:
                if 'size' in log['step']['known_paths'][path]:
                    label += "\\n%s" % misc.bytes_to_str(
                        log['step']['known_paths'][path]['size'])
        pipe_hash['nodes'][misc.str_to_sha1(path)] = {
            'label': label,
            'fillcolor': color
        }
            
    for proc_info in copy.deepcopy(log['pipeline_log']['processes']):
        pid = proc_info['pid']
        # Set name and label variable
        try:
            name = proc_info['name']
            label = "%s" % (proc_info['name'])
        except KeyError:
            name = '(unknown)'
            label = "PID %d" % pid
            
        try:
            # Add file nodes for every file in hints
            for path in proc_info['hints']['writes']:
                add_file_node(path)
                pipe_hash['edges'][(pid_hash(pid), file_hash(path))] = dict()
        except KeyError:
            pass

        try:
            # Add all the info for each process to pipe_hash
            stripped_args = []
            is_output_file = False
            for arg in proc_info['args']:
                # Try to investigate how fifos are integrated in data stream
                # Hier muss die Entscheidung rein ob eine Datei für Input oder
                # Output genutzt wird
                io_type = None
                if name == 'cat':
                    is_output_file = False
                elif name == 'dd' and arg.startswith('of='):
                    is_output_file = True
                elif name == 'dd' and arg.startswith('if='):
                    is_output_file = False
                for known_path in log['step']['known_paths'].keys():
                    # Check if arg contains a known path ...
                    if known_path in arg:
                        # ... if so add this file to the graph 
                        add_file_node(known_path)
                        # Is the process able to in-/output files?
                        if name in ['cat', 'dd']:
                            if is_output_file:
                                io_type = 'output'
                            else:
                                io_type = 'input'
                        else:
                            # we can't know whether the fifo is for input or
                            # output, first look at the hints, then use the
                            # designation (if any was given)
                            if 'reads' in proc_info['hints'] and \
                               arg in proc_info['hints']['reads']:
                                io_type = 'input'
                            if 'writes' in proc_info['hints'] and \
                               arg in proc_info['hints']['writes']:
                                io_type = 'output'
                            if io_type is None:
                                io_type = log['step']['known_paths'][known_path]\
                                            ['designation']
                        if io_type == 'input':
                            # add edge from file to proc
                            pipe_hash['edges']\
                                [(file_hash(known_path),pid_hash(pid))] = dict()
                        elif io_type == 'output':
                            # add edge from proc to file
                            pipe_hash['edges']\
                                [(pid_hash(pid), file_hash(known_path))] = dict()

                        basename = os.path.basename(known_path)
                        #if log['step']['known_paths'][known_path]['type'] != \
                        #   'step_file':
                        arg = arg.replace(known_path, basename)
#                            break
#                        else:
#                            arg = basename
#                            break
                    else:
#                        if arg[0:4] != '/dev':
#                            arg = os.path.basename(arg)
                        if (len(arg) > 16) and re.match('^[A-Z]+$', arg):
                            arg = "%s[...]" % arg[:16]
                stripped_args.append(arg.replace('\t', '\\t').replace(
                    '\\', '\\\\'))
                
            tw = textwrap.TextWrapper(
                width = 50, 
                break_long_words = False, 
                break_on_hyphens = False)
            label = "%s" % ("\\n".join(tw.wrap(' '.join(stripped_args))))

        # If any key wasn't around let's go on
        except KeyError:
            pass

        # add proc
        something_went_wrong = False
        if 'signal' in proc_info:
            something_went_wrong = True
        elif 'exit_code' in proc_info:
            if proc_info['exit_code'] != 0:
                something_went_wrong = True
        else:
            something_went_wrong = True
        color = "#fce94f"
        if something_went_wrong:
            if not pid in log['pipeline_log']['ok_to_fail']:
                color = "#d5291a"
            if 'signal' in proc_info:
                label = "%s\\n(received %s%s)" % (
                    label, 
                    'friendly ' \
                    if pid in log['pipeline_log']['ok_to_fail'] else '',
                    proc_info['signal_name'] if 'signal_name' in \
                    proc_info else 'signal %d' % proc_info['signal'])
            elif 'exit_code' in proc_info:
                if proc_info['exit_code'] != 0:
                    label = "%s\\n(failed with exit code %d)" % (
                        label, proc_info['exit_code'])
            else:
                label = "%s\\n(no exit code)" % label
                    
        if 'max' in log['pipeline_log']['process_watcher']:
            if pid in log['pipeline_log']['process_watcher']['max']:
                label += "\\n%1.1f%% CPU, %s RAM (%1.1f%%)" % (
                    log['pipeline_log']['process_watcher']['max'][pid]\
                    ['cpu_percent'],
                    misc.bytes_to_str(
                        log['pipeline_log']['process_watcher']['max'][pid]\
                        ['rss']),
                    log['pipeline_log']['process_watcher']['max'][pid]\
                    ['memory_percent'])
                
        pipe_hash['nodes'][pid_hash(pid)] = {
            'label': label,
            'fillcolor': color
        }
            
        for which in ['stdout', 'stderr']:
            key = "%s_copy" % which
            if key in proc_info:
                if ('exit_code' in proc_info[key]) and \
                   (proc_info[key]['exit_code'] == 0) and \
                   ('length' in proc_info[key]) and \
                   (proc_info[key]['length'] == 0) and \
                   (not 'sink_full_path' in proc_info[key]):
                    # skip this stdout/stderr box if it leads to nothing
                    continue
                size_label = '(empty)'
                if ('length' in proc_info[key]) and \
                   (proc_info[key]['length'] > 0):
                    speed = float(proc_info[key]['length']) / (
                        proc_info[key]['end_time'] - 
                        proc_info[key]['start_time']).total_seconds()
                    speed_label = "%s/s" % misc.bytes_to_str(speed)
                    size_label = "%s / %s lines (%s)" % (
                        misc.bytes_to_str(proc_info[key]['length']),
                        "{:,}".format(proc_info[key]['lines']),
                        speed_label)
                label = "%s\\n%s" % (which, size_label)
                
                something_went_wrong = False
                if 'signal' in proc_info[key]:
                    something_went_wrong = True
                elif 'exit_code' in proc_info[key]:
                    if proc_info[key]['exit_code'] != 0:
                        something_went_wrong = True
                else:
                    something_went_wrong = True
                color = "#fdf3a7"
                if something_went_wrong:
                    if not pid in log['pipeline_log']['ok_to_fail']:
                        color = "#d5291a"
                    if 'signal' in proc_info[key]:
                        label = "%s\\n(received %s%s)" % (
                            label, 
                            "friendly " if pid in \
                            log['pipeline_log']['ok_to_fail'] else '',
                            proc_info[key]['signal_name'] if 'signal_name'\
                            in proc_info[key] else 'signal %d' % 
                            proc_info[key]['signal'])
                    elif 'exit_code' in proc_info[key]:
                        if proc_info[key]['exit_code'] != 0:
                            label = "%s\\n(failed with exit code %d)" % (
                                label, proc_info[key]['exit_code'])
                    else:
                        label = "%s\\n(no exit code)" % label
                        
                                
                # add proc_which
                pipe_hash['nodes'][pid_hash(pid, which)] = {
                    'label': label,
                    'fillcolor': color
                }
                if 'sink_full_path' in proc_info[key]:
                    path = proc_info[key]['sink_full_path']
                    add_file_node(path)

    for proc_info in copy.deepcopy(log['pipeline_log']['processes']):
        pid = proc_info['pid']
        if 'use_stdin_of' in proc_info:
            other_pid = proc_info['use_stdin_of']
            pipe_hash['edges'][(pid_hash(other_pid, 'stdout'), pid_hash(pid))] \
                = dict()
        for which in ['stdout', 'stderr']:
            key = "%s_copy" % which
            if key in proc_info:
                other_pid = proc_info[key]['pid']
                pipe_hash['edges'][(pid_hash(pid), pid_hash(pid, which))] = dict()
                if 'sink_full_path' in proc_info[key]:
                    pipe_hash['edges'][(
                        pid_hash(pid, which),
                        file_hash(proc_info[key]['sink_full_path']))] = dict()

    # define nodes which go into subgraph
    step_file_nodes = dict()
    for path, path_info in log['step']['known_paths'].items():
        if path_info['type'] == 'step_file':
            step_file_nodes[file_hash(path)] = path_info['designation']

    task_name = "%s/%s" % (log['step']['name'], log['run']['run_id'])
    cluster_hash = misc.str_to_sha1(task_name)
    pipe_hash['clusters'][cluster_hash] = dict()
    pipe_hash['clusters'][cluster_hash]['task_name'] = task_name
    pipe_hash['clusters'][cluster_hash]['group'] = list()
    for node in pipe_hash['nodes'].keys():
        if not node in step_file_nodes:
            pipe_hash['clusters'][cluster_hash]['group'].append(node)
                
    start_time = log['start_time']
    end_time = log['end_time']
    duration = end_time - start_time

    text = "Task: %s\\lHost: %s\\lDuration: %s\\l" % (
        task_name, socket.gethostname(),
        misc.duration_to_str(duration, long = True)
    )
    pipe_hash['graph_labels'][task_name] = text
    if 'max' in log['pipeline_log']['process_watcher']:
        text = "CPU: %1.1f%%, %d CORES_Requested , RAM: %s (%1.1f%%)\\l" % (
            log['pipeline_log']['process_watcher']['max']['sum']['cpu_percent'],
            log['step']['cores'],
            misc.bytes_to_str(log['pipeline_log']['process_watcher']['max']\
                              ['sum']['rss']), 
            log['pipeline_log']['process_watcher']['max']['sum']['memory_percent'])
        pipe_hash['graph_labels'][task_name] += text
    if 'signal' in log:
        pipe_hash['graph_labels'][task_name] += "Caught signal: %s\\l" % (
            process_pool.ProcessPool.SIGNAL_NAMES[log['signal']])
    pipe_hash['graph_labels'][task_name] += "\\l"
    return pipe_hash
