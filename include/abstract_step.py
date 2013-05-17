import sys
sys.path.append('./include/steps')
sys.path.append('./include/sources')
import copy
import datetime
import fscache
import hashlib
import inspect
import json
import misc
import os
import re
import random
import signal
import socket
import string
import StringIO
import subprocess
import tempfile
import textwrap
import time
import traceback
import yaml


class AbstractStep(object):
    
    fsc = fscache.FSCache()
    
    PING_TIMEOUT = 300
    PING_RENEW = 30
    
    def __init__(self, pipeline):
        
        self._pipeline = pipeline
        
        self.dependencies = []
        '''
        All steps this step depends on.
        '''
        
        self.options = {}
        '''
        Options as specified in the configuration.
        '''
        
        self._step_name = self.__module__
        '''
        By default, this is the name of the module. Can be overridden 
        to allow for multiple steps of the same kind.
        '''
        
        self._run_info = None
        '''
        Cached run information. ``setup_runs`` is only called once, the
        post-processed results are stored in here.
        '''
        
        self._temp_directory = None
        '''
        The temporary output directory the step is using. Only set when
        the step is being run.
        '''
        
        self._pipeline_log = dict()
        
        self._cores = 1
        self._connections = set()
        self._connection_restrictions = {}
        self._tools = dict()
        
        self.needs_parents = False
        
        self.known_paths = dict()
        
        self.children_step_names = set()
        
        self.finalized = False
        
    def finalize(self):
        if self.finalized:
            return
        
        # find out which steps are in our family tree
        self.ancestors = set()
        for parent_step in self.dependencies:
            parent_step.finalize()
            self.ancestors.add(str(parent_step))
            for grand_parent_step in parent_step.ancestors:
                self.ancestors.add(grand_parent_step)
            
        self.finalized = True
        
    def _reset(self):
        self.known_paths = dict()
        self._pipeline_log = dict()

    def set_name(self, step_name):
        self._step_name = step_name

    def set_options(self, options):
        self.options = options

    def add_dependency(self, parent):
        if not isinstance(parent, AbstractStep):
            raise StandardError("Error: parent argument must be an AbstractStep.")
        if parent == self:
            raise StandardError("Cannot add a node as its own dependency.")
        self.dependencies.append(parent)
        parent.children_step_names.add(str(self))
        
    def get_input_run_info(self):
        '''
        Return a dict with run info for each parent.
        '''
        input_run_info = dict()
        for parent in self.dependencies:
            input_run_info[parent.get_step_name()] = copy.deepcopy(parent.get_run_info())
        return input_run_info

    def setup_runs(self, input_run_info, connection_info = {}):
        '''
        Raise NotImplementedError because every subclass must override this method.
        '''
        raise NotImplementedError()

    def execute(self, run_id, run_info):
        '''
        Raise NotImplementedError because every subclass must override this method.
        '''
        raise NotImplementedError()

    def get_run_info(self):
        # create run info if it doesn't exist yet
        if not self._run_info:
            # if _BREAK: true is specified in the configuration,
            # return no runs and thus cut off further processing
            if '_BREAK' in self.options and self.options['_BREAK']:
                return dict()
                
            # create input run info and simplify it a bit for setup_runs
            input_run_info = copy.deepcopy(self.get_input_run_info())
            full_paths = dict()

            # strip directories from file names, strip input files
            for step_name in input_run_info.keys():
                for run_id, run_info in input_run_info[step_name].items():
                    for tag in run_info['output_files'].keys():
                        for path in run_info['output_files'][tag].keys():
                            basename = os.path.basename(path)
                            if basename in full_paths:
                                raise StandardError("There are multiple input filenames with the same basename.")
                            full_paths[basename] = path
                            run_info['output_files'][tag][basename] = run_info['output_files'][tag][path]
                            run_info['output_files'][tag].pop(path)

            connection_info = {}
            for c in self._connections:
                if c[0:3] != 'in/':
                    continue
                connection_info[c] = self.get_input_run_info_for_connection(c)
                
            self._run_info = self.setup_runs(input_run_info, connection_info)
            
            # verify run_ids and connection keys
            for run_id, run_info in self._run_info.items():
                if '/' in run_id:
                    raise StandardError("Run IDs must not contain a slash ('/'): %s "
                        "returns a run called %s." % (self, run_id))
                for tag in run_info['output_files'].keys():
                    if not 'out/' + tag in self._connections:
                        raise StandardError("Invalid output_file tag '%s' in %s. "
                            "You might want to add self.add_connection('out/%s')"
                            "to the constructor of %s." 
                            % (tag, str(self), tag, self.__class__))
            
            for run_id, run_info in self._run_info.items():
                for tag in run_info['output_files'].keys():
                    for path in run_info['output_files'][tag].keys():
                        full_paths[path] = os.path.join(self.get_output_directory(), path)

            self._run_info = fix_dict(self._run_info, fix_func_dict_subst, full_paths)
                        
            # define file dependencies
            for run_id in self._run_info.keys():
                for tag in self._run_info[run_id]['output_files'].keys():
                    for output_path, input_paths in self._run_info[run_id]['output_files'][tag].items():
                        self._pipeline.add_file_dependencies(output_path, input_paths)
                        self._pipeline.add_output_file_created_by(output_path, str(self), run_id)
                        
        # now that the _run_info exists, it remains constant, just return it
        return self._run_info

    def get_run_ids(self):
        return self.get_run_info().keys()

    def get_options_hashtag(self):
        options_without_dash_prefix = dict()
        for k, v in self.options.items():
            if k[0] != '_':
                options_without_dash_prefix[k] = v
        return hashlib.sha1(json.dumps(options_without_dash_prefix, sort_keys=True)).hexdigest()[0:4]

    def get_step_name(self):
        return self._step_name

    def get_output_directory(self):
        return os.path.join(self._pipeline.config['destination_path'], 
            '%s-%s' % (self.get_step_name(), self.get_options_hashtag()))

    def get_temp_output_directory(self):
        while True:
            token = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(8))
            path = os.path.join(self._pipeline.config['destination_path'], 'temp', 'temp-%s-%s' % (str(self), token))
            if not os.path.exists(path):
                return path

    def get_run_state(self, run_id):

        def path_up_to_date(outpath, inpaths):
            if not AbstractStep.fsc.exists(outpath):
                return False
            for inpath in inpaths:
                if not AbstractStep.fsc.exists(inpath):
                    return False
                if AbstractStep.fsc.getmtime(inpath) > AbstractStep.fsc.getmtime(outpath):
                    return False
            return True
            
        def up_to_dateness_level(path, level = 0):
            #print("up_to_dateness_level(path = %s, level = %d)" % (os.path.basename(path), level))
            result = level
            dep_paths = self._pipeline.file_dependencies[path]
            if not path_up_to_date(path, dep_paths):
                result = level + 1
            for dep_path in dep_paths:
                recursive_result = up_to_dateness_level(dep_path, level + 1)
                if recursive_result > level + 1:
                    result = max(result, recursive_result)
            return result

        '''
        - finished: all output files exist AND up to date (recursively)
        - ready: NOT all output files exist AND all input files exist AND up to date (recursively)
        - waiting: otherwise
        - if it's ready, it might be executing or queued -> check execute and queue ping
        - if it's waiting, it might be queued -> check queue ping
        
        the ping works like this (this example is for execute, same goes for queued):
          - there's a ping file for every task ( = step + run)
          - it contains information about when how where the job was started etc.
          - its timestamp gets renewed every 30 seconds (touch)
          - as soon as the job has finished, the execute ping file is removed,
            this should also work if the job crashes (however, it cannot work if
            the controlling script receives SIGKILL
          - if its timestamp is no more than 5 minutes old, it is regarded as
            currently executing
          - otherwise, a warning is printed because the ping file is probably stale
            (no automatic cleanup is performed, manual intervention is necessary)
          - warning: this requires all involved systems or the file system to be
            time-synchronized
        '''
        
        run_info = self.get_run_info()
        max_level = 0
        for tag in run_info[run_id]['output_files'].keys():
            for output_file in run_info[run_id]['output_files'][tag].keys():
                max_level = max(max_level, up_to_dateness_level(output_file))

        if max_level == 0:
            return self._pipeline.states.FINISHED
        elif max_level == 1:
            executing_ping_path = self.get_executing_ping_path_for_run_id(run_id)
            if AbstractStep.fsc.exists(executing_ping_path):
                return self._pipeline.states.EXECUTING
            if AbstractStep.fsc.exists(self.get_queued_ping_path_for_run_id(run_id)):
                return self._pipeline.states.QUEUED
            return self._pipeline.states.READY
        else:
            if AbstractStep.fsc.exists(self.get_queued_ping_path_for_run_id(run_id)):
                return self._pipeline.states.QUEUED
            return self._pipeline.states.WAITING

    def run(self, run_id):
        # create the output directory if it doesn't exist yet
        if not os.path.isdir(self.get_output_directory()):
            os.makedirs(self.get_output_directory())
            
        # now write the run ping file
        executing_ping_path = self.get_executing_ping_path_for_run_id(run_id)
        
        if os.path.exists(executing_ping_path):
            raise StandardError("%s/%s seems to be already running, exiting..." % (self, run_id))
        
        queued_ping_path = self.get_queued_ping_path_for_run_id(run_id)
        
        # create a temporary directory for the output files
        temp_directory = self.get_temp_output_directory()
        self._temp_directory = temp_directory
        os.makedirs(temp_directory)

        # call execute() but pass output file paths with the temporary directory
        temp_run_info = copy.deepcopy(self.get_run_info()[run_id])

        # prepare self.known_paths
        self.known_paths = dict()
        for tag, tag_info in temp_run_info['output_files'].items():
            for output_path, input_paths in tag_info.items():
                # add the real output path
                self.known_paths[output_path] = {'type': 'output', 'designation': 'output', 'label': os.path.basename(output_path), 'type': 'step_file'}
                # ...and also add the temporary output path
                self.known_paths[os.path.join(temp_directory, os.path.basename(output_path))] = {'type': 'output', 'designation': 'output', 'label': "%s\\n(%s)" % (os.path.basename(output_path), tag), 'type': 'step_file', 'real_path': output_path}
                for input_path in input_paths:
                    self.known_paths[input_path] = {'type': 'input', 'designation': 'input', 'label': os.path.basename(input_path), 'type': 'step_file'}

        temp_paths = {}
        for tag in temp_run_info['output_files'].keys():
            for out_path, in_paths in temp_run_info['output_files'][tag].items():
                temp_paths[out_path] = os.path.join(temp_directory, os.path.basename(out_path))

        temp_run_info = fix_dict(temp_run_info, fix_func_dict_subst, temp_paths)
        
        # now write the run ping file
        executing_ping_info = dict()
        executing_ping_info['start_time'] = datetime.datetime.now()
        executing_ping_info['host'] = socket.gethostname()
        executing_ping_info['pid'] = os.getpid()
        executing_ping_info['cwd'] = os.getcwd()
        with open(executing_ping_path, 'w') as f:
            f.write(yaml.dump(executing_ping_info, default_flow_style = False))
            
        executing_ping_pid = os.fork()
        if executing_ping_pid == 0:
            while True:
                time.sleep(AbstractStep.PING_RENEW)
                if os.path.exists(executing_ping_path):
                    os.utime(executing_ping_path, None)
            os._exit(0)
            
        self.start_time = datetime.datetime.now()
        self._pipeline.notify("[INFO] starting %s/%s on %s" % (str(self), run_id, socket.gethostname()))
        try:
            self.execute(run_id, temp_run_info)
        except:
            self.end_time = datetime.datetime.now()
            annotation_path, annotation_str = self.write_annotation(run_id, self._temp_directory)
            message = "[BAD] %s/%s failed on %s\n\nHere are the details:\n%s" % (str(self), run_id, socket.gethostname(), annotation_str)
            attachment = None
            if os.path.exists(annotation_path + '.png'):
                attachment = dict()
                attachment['name'] = 'details.png'
                attachment['data'] = open(annotation_path + '.png').read()
            self._pipeline.notify(message, attachment)
            raise
        finally:
            try:
                os.kill(executing_ping_pid, signal.SIGTERM)
                os.waitpid(executing_ping_pid, 0)
            except OSError:
                # if the ping process was already killed, it's gone anyway
                pass
            # remove the run ping file
            if os.path.exists(executing_ping_path):
                try:
                    os.unlink(executing_ping_path)
                except OSError:
                    pass
            # remove the queued ping file
            if os.path.exists(queued_ping_path):
                try:
                    os.unlink(queued_ping_path)
                except OSError:
                    pass
                
        # TODO: Clean this up. Re-think exceptions and task state transisitions.
        
        self.end_time = datetime.datetime.now()
        
        # if we're here, we can assume the step has finished successfully
        # now rename the output files (move from temp directory to
        # destination directory)
        for tag in temp_run_info['output_files'].keys():
            for out_path in temp_run_info['output_files'][tag].keys():
                destination_path = os.path.join(self.get_output_directory(), os.path.basename(out_path))
                # TODO: if the destination path already exists, this will overwrite the file.
                if os.path.exists(out_path):
                    os.rename(out_path, destination_path)
                else:
                    annotation_path, annotation_str = self.write_annotation(run_id, self._temp_directory)
                    raise StandardError("The step failed to produce an output file it announced: %s\n\nHere are the details:\n%s" % (os.path.basename(out_path), annotation_str))

        for path, path_info in self.known_paths.items():
            if os.path.exists(path):
                self.known_paths[path]['size'] = os.path.getsize(path)

        self._temp_directory = None

        annotation_path, annotation_str = self.write_annotation(run_id, self.get_output_directory())

        # create a symbolic link to the annotation for every output file
        for tag in temp_run_info['output_files'].keys():
            for out_path in temp_run_info['output_files'][tag].keys():
                destination_path = os.path.join(self.get_output_directory(), '.' + os.path.basename(out_path) + '.annotation.yaml')
                # overwrite the symbolic link if it already exists
                if os.path.exists(destination_path):
                    os.unlink(destination_path)
                oldwd = os.getcwd()
                os.chdir(os.path.dirname(destination_path))
                os.symlink(os.path.basename(annotation_path), os.path.basename(destination_path))
                os.chdir(oldwd)

        # finally, remove the temporary directory if it's empty
        try:
            os.rmdir(temp_directory)
        except OSError:
            pass
        
        # step has completed successfully, now determine how many jobs are still left
        # but first invalidate the FS cache because things have changed by now...
        AbstractStep.fsc = fscache.FSCache()
        
        remaining_task_info = self.get_run_info_str()

        message = "[OK] %s/%s successfully finished on %s\n" % (str(self), run_id, socket.gethostname())
        message += str(self) + ': ' + remaining_task_info + "\n"
        attachment = None
        if os.path.exists(annotation_path + '.png'):
            attachment = dict()
            attachment['name'] = 'details.png'
            attachment['data'] = open(annotation_path + '.png').read()
        self._pipeline.notify(message, attachment)

        self._reset()

    def tool(self, key):
        '''
        Return full path to a configured tool.
        '''
        return copy.deepcopy(self._tools[key])
    
    def get_run_info_str(self):
        count = {}
        for _ in self.get_run_ids():
            state = self.get_run_state(_)
            if not state in count:
                count[state] = 0
            count[state] += 1
        return ', '.join(["%d %s" % (count[_], _.lower()) for _ in self._pipeline.states.order if _ in count])
    
    def append_pipeline_log(self, log):
        if len(self._pipeline_log) == 0:
            self._pipeline_log = copy.deepcopy(log)
        else:
            for k in log.keys():
                if k == 'process_watcher':
                    for k2 in log[k].keys():
                        if k2 == 'max':
                            for _ in log[k][k2].keys():
                                if _ == 'sum':
                                    for k3 in self._pipeline_log[k][k2][_].keys():
                                        self._pipeline_log[k][k2][_][k3] = max(self._pipeline_log[k][k2][_][k3], log[k][k2][_][k3])
                                else:
                                    self._pipeline_log[k][k2][_] = copy.deepcopy(log[k][k2][_])
                        else:
                            self._pipeline_log[k][k2].update(log[k][k2])
                            
                else:
                    if log[k].__class__ == list:
                        self._pipeline_log[k].extend(log[k])
                    else:
                        self._pipeline_log[k].update(log[k])
    
    def write_annotation(self, run_id, path):
        # now write the annotation
        log = {}
        log['step'] = {}
        log['step']['options'] = self.options
        log['step']['name'] = self.get_step_name()
        log['step']['known_paths'] = self.known_paths
        log['run'] = {}
        log['run']['run_info'] = self.get_run_info()[run_id]
        log['run']['run_id'] = run_id
        log['config'] = self._pipeline.config
        log['git_hash_tag'] = self._pipeline.git_hash_tag
        log['tool_versions'] = {}
        for tool in self._tools.keys():
            log['tool_versions'][tool] = self._pipeline.tool_versions[tool]
        log['pipeline_log'] = self._pipeline_log
        log['start_time'] = self.start_time
        log['end_time'] = self.end_time
        if self._pipeline.git_dirty_diff:
            log['git_dirty_diff'] = self._pipeline.git_dirty_diff

        annotation_yaml = yaml.dump(log, default_flow_style = False)
        annotation_path = os.path.join(path, ".%s-annotation-%s.yaml" % (run_id, hashlib.sha1(annotation_yaml).hexdigest()[0:8]))
        
        # overwrite the annotation if it already exists
        with open(annotation_path, 'w') as f:
            f.write(annotation_yaml)
            
        try:
            gv = self.render_pipeline([log])
            dot = subprocess.Popen(['dot', '-Tsvg'], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
            dot.stdin.write(gv)
            dot.stdin.close()
            svg = dot.stdout.read()
            with open(annotation_path + '.svg', 'w') as f:
                f.write(svg)
                
            dot = subprocess.Popen(['dot', '-Tpng'], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
            dot.stdin.write(gv)
            dot.stdin.close()
            png = dot.stdout.read()
            with open(annotation_path + '.png', 'w') as f:
                f.write(png)
        except:
            # rendering the pipeline graph is not _that_ important, after all
            # we can still try to render it later from the annotation file
            pass
        
        return annotation_path, annotation_yaml
    
    def get_temporary_path(self, prefix = '', designation = None):
        '''
        Returns a temporary path with the prefix specified. 
        The returned path will be in the temporary directory of the step 
        and will not exist yet.
        '''
        if not self._temp_directory:
            raise StandardError("Temporary directory not set, you cannot call get_temporary_path from setup_runs.")

        _, _path = tempfile.mkstemp('', prefix, self._temp_directory)
        os.close(_)
        os.unlink(_path)
        
        self.known_paths[_path] = {'label': prefix, 'designation': designation, 'type': 'file'}

        return _path        

    def get_temporary_fifo(self, prefix = '', designation = None):
        '''
        Create a temporary FIFO and return its path.
        '''
        path = self.get_temporary_path(prefix, designation)
        os.mkfifo(path)
        self.known_paths[path]['type'] = 'fifo'
        return path

    def __str__(self):
        return self._step_name
    
    @classmethod
    def render_pipeline(cls, logs):
        hash = {'nodes': {}, 'edges': {}, 'clusters': {}, 'graph_labels': {}}
        for log in logs:
            temp = cls.render_pipeline_hash(log)
            for _ in ['nodes', 'edges', 'clusters', 'graph_labels']:
                hash[_].update(temp[_])

        f = StringIO.StringIO()
        f.write("digraph {\n")
        f.write("    rankdir = TB;\n")
        f.write("    splines = true;\n")
        f.write("    graph [fontname = Helvetica, fontsize = 12, size = \"14, 11\", nodesep = 0.2, ranksep = 0.3, labelloc = t, labeljust = l];\n")
        f.write("    node [fontname = Helvetica, fontsize = 12, shape = rect, style = filled];\n")
        f.write("    edge [fontname = Helvetica, fontsize = 12];\n")
        f.write("\n")
        
        f.write("    // nodes\n")
        f.write("\n")
        for node_key, node_info in hash['nodes'].items():
            f.write("    _%s" % node_key)
            if len(node_info) > 0:
                f.write(" [%s]" % ', '.join(['%s = "%s"' % (k, node_info[k]) for k in node_info.keys()]))
            f.write(";\n")
            
        f.write("\n")
        
        f.write("    // edges\n")
        f.write("\n")
        for edge_pair in hash['edges'].keys():
            if edge_pair[0] in hash['nodes'] and edge_pair[1] in hash['nodes']:
                f.write("    _%s -> _%s;\n" % (edge_pair[0], edge_pair[1]))
        
        f.write("\n")
        
        '''
        f.write("    // clusters\n")
        f.write("\n")
        for cluster_hash, cluster_info in hash['clusters'].items():
            f.write("    subgraph cluster_%s {\n" % cluster_hash)
            for node in cluster_info['group']:
                f.write("        _%s;\n" % node)
                
            f.write("        label = \"%s\";\n" % cluster_info['task_name'])
            f.write("        graph [style = dashed];\n")
            f.write("    }\n")
        '''
        
        if len(hash['graph_labels']) == 1:
            f.write("    graph [label=\"%s\"];\n" % hash['graph_labels'].values()[0])
        f.write("}\n")
        
        result = f.getvalue()
        f.close()
        return result
        
    @classmethod
    def render_pipeline_hash(cls, log):
        
        def pid_hash(pid, suffix = ''):
            hashtag = "%s/%s/%d/%s" % (log['step']['name'], log['run']['run_id'], pid, suffix)
            return hashlib.sha1(hashtag).hexdigest()
        
        def file_hash(path):
            if 'real_path' in log['step']['known_paths'][path]:
                path = log['step']['known_paths'][path]['real_path']
            return misc.str_to_sha1(path)
        
        #print(yaml.dump(self.known_paths, default_flow_style = False))
        
        hash = dict()
        hash['nodes'] = dict()
        hash['edges'] = dict()
        hash['clusters'] = dict()
        hash['graph_labels'] = dict()
        
        def add_file_node(path):
            if 'real_path' in log['step']['known_paths'][path]:
                path = log['step']['known_paths'][path]['real_path']
            label = log['step']['known_paths'][path]['label']
            color = '#ffffff'
            if log['step']['known_paths'][path]['type'] == 'fifo':
                color = '#c4f099'
            elif log['step']['known_paths'][path]['type'] == 'file':
                color = '#8ae234'
            elif log['step']['known_paths'][path]['type'] == 'step_file':
                color = '#97b7c8'
                if path in log['step']['known_paths']:
                    if 'size' in log['step']['known_paths'][path]:
                        label += "\\n%s" % misc.bytes_to_str(log['step']['known_paths'][path]['size'])
            hash['nodes'][misc.str_to_sha1(path)] = {
                'label': label,
                'fillcolor': color
            }
            
        
        for proc_info in copy.deepcopy(log['pipeline_log']['processes']):
            pid = proc_info['pid']
            label = "PID %d" % pid
            name = '(unknown)'
            if 'name' in proc_info:
                name = proc_info['name']
            label = "%s" % (proc_info['name'])
            if 'writes' in proc_info['hints']:
                for path in proc_info['hints']['writes']:
                    add_file_node(path)
            if 'args' in proc_info:
                stripped_args = []
                for arg in copy.deepcopy(proc_info['args']):
                    if arg in log['step']['known_paths']:
                        add_file_node(arg)
                    if arg in log['step']['known_paths']:
                        if log['step']['known_paths'][arg]['type'] != 'step_file':
                            arg = log['step']['known_paths'][arg]['label']
                        else:
                            arg = os.path.basename(arg)
                    else:
                        if arg[0:4] != '/dev':
                            arg = os.path.basename(arg)
                            if (len(arg) > 16) and re.match('^[A-Z]+$', arg):
                                arg = "%s[...]" % arg[:16]
                    stripped_args.append(arg.replace('\t', '\\t').replace('\\', '\\\\'))
                tw = textwrap.TextWrapper(width = 50, break_long_words = False, break_on_hyphens = False)
                label = "%s" % ("\\n".join(tw.wrap(' '.join(stripped_args))))
            if 'args' in proc_info:
                cat4m_seen_minus_o = False
                for arg in proc_info['args']:
                    fifo_type = None
                    if name == 'cat4m' and arg == '-o':
                        cat4m_seen_minus_o = True
                    if arg in log['step']['known_paths']:
                        add_file_node(arg)
                        if name == 'cat4m':
                            if cat4m_seen_minus_o:
                                fifo_type = 'output'
                            else:
                                fifo_type = 'input'
                        else:
                            # we can't know whether the fifo is for input or output,
                            # firts look at the hints, then use the designation (if any was given)
                            if 'reads' in proc_info['hints'] and arg in proc_info['hints']['reads']:
                                fifo_type = 'input'
                            if 'writes' in proc_info['hints'] and arg in proc_info['hints']['writes']:
                                fifo_type = 'output'
                            if fifo_type is None:
                                fifo_type = log['step']['known_paths'][arg]['designation']
                        if fifo_type == 'input':
                            # add edge from file to proc
                            hash['edges'][(file_hash(arg), pid_hash(pid))] = dict()
                        elif fifo_type == 'output':
                            # add edge from proc to file
                            hash['edges'][(pid_hash(pid), file_hash(arg))] = dict()
            if 'writes' in proc_info['hints']:
                for path in proc_info['hints']['writes']:
                    hash['edges'][(pid_hash(pid), file_hash(path))] = dict()
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
                color = "#d5291a"
                if 'signal' in proc_info:
                    label = "%s\\n(received %s)" % (label, proc_info['signal_name'] if 'signal_name' in proc_info else 'signal %d' % proc_info['signal'])
                elif 'exit_code' in proc_info:
                    if proc_info['exit_code'] != 0:
                        label = "%s\\n(failed with exit code %d)" % (label, proc_info['exit_code'])
                else:
                    label = "%s\\n(exited instantly)" % label
                    
            if 'max' in log['pipeline_log']['process_watcher']:
                if pid in log['pipeline_log']['process_watcher']['max']:
                    label += "\\n%1.1f%% CPU, %s RAM (%1.1f%%)" % (log['pipeline_log']['process_watcher']['max'][pid]['cpu_percent'], misc.bytes_to_str(log['pipeline_log']['process_watcher']['max'][pid]['rss']), log['pipeline_log']['process_watcher']['max'][pid]['memory_percent'])
                
            hash['nodes'][pid_hash(pid)] = {
                'label': label,
                'fillcolor': color
            }
            
            for which in ['stdout', 'stderr']:
                key = "%s_copy" % which
                if key in proc_info:
                    if ('exit_code' in proc_info[key]) and (proc_info[key]['exit_code'] == 0) and ('length' in proc_info[key]) and (proc_info[key]['length'] == 0) and (not 'sink_full_path' in proc_info[key]):
                        # skip this stdout/stderr box if it leads to nothing
                        continue
                    size_label = '(empty)'
                    if ('length' in proc_info[key]) and (proc_info[key]['length'] > 0):
                        speed = float(proc_info[key]['length']) / (proc_info[key]['end_time'] - proc_info[key]['start_time']).total_seconds()
                        speed_label = "%s/s" % misc.bytes_to_str(speed)
                        size_label = "%s / %s lines (%s)" % (misc.bytes_to_str(proc_info[key]['length']), "{:,}".format(proc_info[key]['lines']), speed_label)
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
                        color = "#d5291a"
                        if 'signal' in proc_info[key]:
                            label = "%s\\n(received %s)" % (label, proc_info[key]['signal_name'] if 'signal_name' in proc_info[key] else 'signal %d' % proc_info[key]['signal'])
                        elif 'exit_code' in proc_info[key]:
                            if proc_info[key]['exit_code'] != 0:
                                label = "%s\\n(failed with exit code %d)" % (label, proc_info[key]['exit_code'])
                        else:
                            label = "%s\\n(exited instantly)" % label
                            
                                
                    # add proc_which
                    hash['nodes'][pid_hash(pid, which)] = {
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
                hash['edges'][(pid_hash(other_pid, 'stdout'), pid_hash(pid))] = dict()
            for which in ['stdout', 'stderr']:
                key = "%s_copy" % which
                if key in proc_info:
                    other_pid = proc_info[key]['pid']
                    hash['edges'][(pid_hash(pid), pid_hash(pid, which))] = dict()
                    if 'sink_full_path' in proc_info[key]:
                        hash['edges'][(pid_hash(pid, which), file_hash(proc_info[key]['sink_full_path']))] = dict()

        # define nodes which go into subgraph
        step_file_nodes = dict()
        for path, path_info in log['step']['known_paths'].items():
            if path_info['type'] == 'step_file':
                step_file_nodes[file_hash(path)] = path_info['designation']

        task_name = "%s/%s" % (log['step']['name'], log['run']['run_id'])
        cluster_hash = misc.str_to_sha1(task_name)
        hash['clusters'][cluster_hash] = dict()
        hash['clusters'][cluster_hash]['task_name'] = task_name
        hash['clusters'][cluster_hash]['group'] = list()
        for node in hash['nodes'].keys():
            if not node in step_file_nodes:
                hash['clusters'][cluster_hash]['group'].append(node)
                
        start_time = log['start_time']
        end_time = log['end_time']
        duration = end_time - start_time
        
        if 'max' in log['pipeline_log']['process_watcher']:
            hash['graph_labels'][task_name] = "Task: %s\\lHost: %s, CPU: %1.1f%% , RAM: %s (%1.1f%%)\\lDuration: %s\\l\\l" % (
                task_name, 
                socket.gethostname(),
                log['pipeline_log']['process_watcher']['max']['sum']['cpu_percent'], 
                misc.bytes_to_str(log['pipeline_log']['process_watcher']['max']['sum']['rss']), 
                log['pipeline_log']['process_watcher']['max']['sum']['memory_percent'],
                duration)
        else:
            hash['graph_labels'][task_name] = "Task: %s\\lHost: %s\\lDuration: %s\\l\\l" % (
                task_name, 
                socket.gethostname(),
                duration)
        return hash

    @classmethod
    def get_step_class_for_key(cls, key):
        '''
        Returns a step (or source step) class for a given key which corresponds
        to the name of the module the class is defined in. Pass 'cutadapt' and
        you will get the cutadapt.Cutadapt class which you may then instantiate.
        '''
        
        # look for a subclass of AbstractSourceStep fist
        classes = [_ for _ in inspect.getmembers(__import__(key), inspect.isclass) if AbstractSourceStep in _[1].__bases__]
        if len(classes) > 0:
            if len(classes) != 1:
                raise StandardError("need exactly one subclass of AbstractSourceStep in " + key)
            return classes[0][1]

        # then, look for a subclass of AbstractStep fist
        classes = [_ for _ in inspect.getmembers(__import__(key), inspect.isclass) if AbstractStep in _[1].__bases__]
        classes = [_ for _ in classes if _[1] != AbstractSourceStep]
        if len(classes) != 1:
            print(classes)
            raise StandardError("need exactly one subclass of AbstractStep in " + key)
        return classes[0][1]
    
    def set_cores(self, cores):
        self._cores = cores

    def add_connection(self, connection, constraints = None):
        if connection[0:3] == 'in/':
            self.needs_parents = True
        self._connections.add(connection)
        if constraints is not None:
            self._connection_restrictions[connection] = constraints
        
    def require_tool(self, tool):
        if self._pipeline is not None:
            if not tool in self._pipeline.config['tools']:
                raise StandardError("%s requires %s but it's not declared in the configuration." % (self, tool))
            self._tools[tool] = copy.deepcopy(self._pipeline.config['tools'][tool]['path'])
        else:
            self._tools[tool] = True
        
    def _get_ping_path_for_run_id(self, run_id, key):
        return os.path.join(self.get_output_directory(), '.%s-%s-ping.yaml' % (run_id, key))
    
    def get_executing_ping_path_for_run_id(self, run_id):
        return self._get_ping_path_for_run_id(run_id, 'run')
        
    def get_queued_ping_path_for_run_id(self, run_id):
        return self._get_ping_path_for_run_id(run_id, 'queued')
        
    def get_input_run_info_for_connection(self, in_key):
        if in_key[0:3] != 'in/':
            raise StandardError("in_key does not start with 'in/': %s" % in_key)
        out_key = in_key.replace('in/', 'out/')
        allowed_steps = None
        if '_connect' in self.options:
            if in_key in self.options['_connect']:
                declaration = self.options['_connect'][in_key]
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
        bare_out_key = out_key.replace('out/', '')
        
        result = dict()
        result['counts'] = {
            'total_steps': 0,
            'total_runs': 0,
            'total_files': 0,
            'min_steps_per_run': None,
            'max_steps_per_run': None,
            'min_files_per_step_and_run': None,
            'max_files_per_step_and_run': None,
            'min_files_per_run': None,
            'max_files_per_run': None,
        }
        
        def update_min_max(key, value):
            for mkey in ['min', 'max']:
                key2 = '%s_%s' % (mkey, key)
                if result['counts'][key2] is None:
                    result['counts'][key2] = value
                result['counts'][key2] = (min if mkey == 'min' else max)(result['counts'][key2], value)
            
        result['runs'] = dict()
        for step_name, step_info in self.get_input_run_info().items():
            if allowed_steps is not None:
                if not step_name in allowed_steps:
                    continue
            if out_key in self._pipeline.steps[step_name]._connections:
                result['counts']['total_steps'] += 1
                for run_id, run_info in step_info.items():
                    result['counts']['total_runs'] += 1
                    paths = run_info['output_files'][bare_out_key].keys()
                    result['counts']['total_files'] += len(paths)
                    if not run_id in result['runs']:
                        result['runs'][run_id] = dict()
                    result['runs'][run_id][step_name] = paths

                    steps_per_run = len(result['runs'][run_id])
                    update_min_max('steps_per_run', steps_per_run)

                    files_per_step_and_run = len(result['runs'][run_id][step_name])
                    update_min_max('files_per_step_and_run', files_per_step_and_run)
                    
                    files_per_run = 0
                    for _ in result['runs'][run_id].values():
                        files_per_run += len(_)
                    update_min_max('files_per_run', files_per_run)
                    
        # check constraints, if any
        if in_key in self._connection_restrictions:
            for k, v in self._connection_restrictions[in_key].items():
                if result['counts'][k] != v:
                    raise StandardError("Connection constraint failed: %s/%s/%s should be %d but is %s." % (self, in_key, k, v, str(result['counts'][k])))

        return result

    
class AbstractSourceStep(AbstractStep):
    '''
    A subclass all source steps inherit from and which distinguishes source
    steps from all real processing steps because they do not yield any tasks, 
    because their "output files" are in fact files which are already there.
    
    Note that the name might be a bit misleading because this class only
    applies to source steps which 'serve' existing files. A step which has 
    no input but produces input data for other steps and actually has to do 
    something for it, on the other hand, would be a normal AbstractStep
    subclass because it produces tasks.
    '''

    def __init__(self, pipeline):
        super(AbstractSourceStep, self).__init__(pipeline)

def fix_dict(data, fix_func, *args):
    if data.__class__ == list:
        return [fix_dict(_, fix_func, *args) for _ in data]
    elif data.__class__ == dict:
        result = {}
        for k, v in data.items():
            result[fix_dict(k, fix_func, *args)] = fix_dict(v, fix_func, *args)
        return result
    else:
        return fix_func(data, *args)

def fix_func_dict_subst(v, full_paths):
    if v in full_paths:
        return full_paths[v]
    return v
