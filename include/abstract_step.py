import sys
sys.path.append('./include/steps')
sys.path.append('./include/sources')
import copy
import datetime
import fscache
import hashlib
import inspect
import json
import os
import random
import socket
import string
import subprocess
import tempfile
import traceback
import unix_pipeline
import yaml


class AbstractStep(object):
    
    fsc = fscache.FSCache()
    
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
        
        self._cores = 1
        self._connections = []
        self._tools = dict()
        
        self.needs_parents = False

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
        
    def get_input_run_info(self):
        '''
        Return a dict with run info for each parent.
        '''
        input_run_info = dict()
        for parent in self.dependencies:
            input_run_info[parent.get_step_name()] = copy.deepcopy(parent.get_run_info())
        return input_run_info

    def setup_runs(self, input_run_info):
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

            self._run_info = self.setup_runs(input_run_info)
            
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
                        
            # fill _file_dependencies
            for run_id in self._run_info.keys():
                for tag in self._run_info[run_id]['output_files'].keys():
                    for output_path, input_paths in self._run_info[run_id]['output_files'][tag].items():
                        self._pipeline.add_file_dependencies(output_path, input_paths)
                        
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
        finished: all output files exist AND up to date (recursively)
        ready: NOT all output files exist AND all input files exist AND up to date (recursively)
        waiting: otherwise
        '''
        
        run_info = self.get_run_info()
        max_level = 0
        for tag in run_info[run_id]['output_files'].keys():
            for output_file in run_info[run_id]['output_files'][tag].keys():
                max_level = max(max_level, up_to_dateness_level(output_file))

        if max_level == 0:
            return self._pipeline.states.FINISHED
        elif max_level == 1:
            return self._pipeline.states.READY
        else:
            return self._pipeline.states.WAITING

    def run(self, run_id):
        # create the output directory if it doesn't exist yet
        if not os.path.isdir(self.get_output_directory()):
            os.makedirs(self.get_output_directory())
            
        # also create a temporary directory for the output file
        temp_directory = self.get_temp_output_directory()
        self._temp_directory = temp_directory
        os.makedirs(temp_directory)

        # call execute() but pass output file paths with the temporary directory
        temp_run_info = copy.deepcopy(self.get_run_info()[run_id])
        temp_paths = {}
        for tag in temp_run_info['output_files'].keys():
            for out_path, in_paths in temp_run_info['output_files'][tag].items():
                temp_paths[out_path] = os.path.join(temp_directory, os.path.basename(out_path))

        temp_run_info = fix_dict(temp_run_info, fix_func_dict_subst, temp_paths)

        start_time = datetime.datetime.now()
        self._pipeline.notify("[INFO] starting %s/%s on %s" % (str(self), run_id, socket.gethostname()))
        try:
            self.execute(run_id, temp_run_info)
        except Exception as e:
            self._pipeline.notify("[BAD] %s/%s failed.\n\nHere are the details:\n%s" % (str(self), run_id, yaml.dump(unix_pipeline.get_log(), default_flow_style = False)))
            raise

        # if we're here, we can assume the step has finished successfully
        # now rename the output files (move from temp directory to
        # destination directory)
        for tag in temp_run_info['output_files'].keys():
            for out_path in temp_run_info['output_files'][tag].keys():
                destination_path = os.path.join(self.get_output_directory(), os.path.basename(out_path))
                # TODO: if the destination path already exists, this will overwrite the file.
                if not os.path.exists(out_path):
                    raise StandardError("The step failed to produce an output file it announced: " + os.path.basename(out_path))
                os.rename(out_path, destination_path)

        self._temp_directory = None
        
        # step has completed successfully, now determine how many jobs are still left
        # but first invalidate the FS cache because things have changed by now...
        AbstractStep.fsc = fscache.FSCache()
        
        count = {}
        for _ in self.get_run_ids():
            state = self.get_run_state(_)
            if not state in count:
                count[state] = 0
            count[state] += 1
        remaining_task_info = ', '.join([str(count[_]) + ' ' + _.lower() for _ in sorted(count.keys())])

        end_time = datetime.datetime.now()
        message = "[OK] %s/%s successfully finished.\n" % (str(self), run_id)
        message += str(self) + ': ' + remaining_task_info + "\n"
        self._pipeline.notify(message)

        # now write the annotation
        log = {}
        log['step'] = {}
        log['step']['options'] = self.options
        log['step']['name'] = self.get_step_name()
        log['run'] = {}
        log['run']['run_info'] = self.get_run_info()[run_id]
        log['run']['run_id'] = run_id
        log['config'] = self._pipeline.config
        log['git_hash_tag'] = self._pipeline.git_hash_tag
        log['tool_versions'] = {}
        for tool in self._tools.keys():
            log['tool_versions'][tool] = self._pipeline.tool_versions[tool]
        log['pipeline_log'] = unix_pipeline.get_log()
        log['start_time'] = start_time
        log['end_time'] = end_time

        annotation_path = os.path.join(self.get_output_directory(), '.' + run_id + '-annotation.yaml')
        # overwrite the annotation if it already exists
        with open(annotation_path, 'w') as f:
            f.write(yaml.dump(log, default_flow_style = False))

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
        
        # finally, reset the unix pipeline module
        unix_pipeline.clear()

    def tool(self, key):
        '''
        Return full path to a configured tool.
        '''
        return self._tools[key]
    
    def get_temporary_path(self, prefix, suffix):
        '''
        Returns a temporary path with the prefix and suffix specified. 
        The returned path will be in the temporary directory of the step 
        and will not exist yet.
        '''
        if not self._temp_directory:
            raise StandardError("Temporary directory not set, you cannot call get_temporary_path from setup_runs.")

        _, _path = tempfile.mkstemp(suffix, prefix, self._temp_directory)
        os.close(_)
        os.unlink(_path)

        return _path        

    def __str__(self):
        return self._step_name

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

    def add_connection(self, connection):
        if connection[0:3] == 'in/':
            self.needs_parents = True
        self._connections.append(connection)
        
    def require_tool(self, tool):
        self._tools[tool] = self._pipeline.config['tools'][tool]['path']
    
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
