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
import string
import subprocess
import tempfile
import traceback
import unix_pipeline
import yaml


class AbstractStep(object):
    
    cores = 1
    connections = []
    
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
        
        self._file_dependencies_cumulative = {}
        '''
        _file_dependencies_cumulative is secondary data which gets generated 
        from _run_info and keeps track of each output file's file 
        dependencies (regardless of the run id) and by file dependencies
        I mean ALL file dependencies including those from all of its
        parent steps (yup, that may become a lot of files).
        '''
        
        self._temp_directory = None
        '''
        The temporary output directory the step is using. Only set when
        the step is being run.
        '''

    def set_name(self, step_name):
        self._step_name = step_name

    def set_options(self, options):
        self.options = options

    def add_dependency(self, parent):
        if not isinstance(parent, AbstractStep):
            raise StandardError("Error: parent argument must be an AbstractStep.")
        if parent == self:
            raise StandardError("Cannot add a node as its own dependency.")
        # TODO: Check for cycles.
        self.dependencies.append(parent)

    def get_input_run_info(self):
        '''
        Return a dict with run info for each parent.
        '''
        input_run_info = dict()
        for parent in self.dependencies:
            input_run_info[parent.get_step_id()] = copy.deepcopy(parent.get_run_info())
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
            # create input run info and simplify it a bit for setup_runs
            input_run_info = copy.deepcopy(self.get_input_run_info())
            full_paths = dict()

            # strip directories from file names, strip input files
            for step_name in input_run_info.keys():
                for run_id, run_info in input_run_info[step_name].items():
                    for annotation in run_info['output_files'].keys():
                        for path in run_info['output_files'][annotation].keys():
                            basename = os.path.basename(path)
                            if basename in full_paths:
                                raise StandardError("There are multiple input filenames with the same basename.")
                            full_paths[basename] = path
                            run_info['output_files'][annotation][basename] = run_info['output_files'][annotation][path]
                            run_info['output_files'][annotation].pop(path)

            self._run_info = self.setup_runs(input_run_info)
            
            # verify connection keys
            for run_id, run_info in self._run_info.items():
                for annotation in run_info['output_files'].keys():
                    if not 'out/' + annotation in self.__class__.connections:
                        raise StandardError("Invalid output_file key '%s' in %s. "
                            "Keys must be specified via the 'connections' "
                            "class member (you might want to add "
                            "connections.append('out/%s'))." 
                            % (annotation, str(self), annotation))
            
            for run_id, run_info in self._run_info.items():
                for annotation in run_info['output_files'].keys():
                    for path in run_info['output_files'][annotation].keys():
                        full_paths[path] = os.path.join(self.get_output_directory(), path)

            self._run_info = fix_dict(self._run_info, fix_func_dict_subst, full_paths)
            
            # TODO: Fix this.
            # fill _file_dependencies_cumulative
            for run_id in self._run_info.keys():
                for tag in self._run_info[run_id]['output_files'].keys():
                    for output_path, input_paths in self._run_info[run_id]['output_files'][tag].items():
                        if output_path in self._file_dependencies_cumulative:
                            raise StandardError("Error: Multiple run IDs want to create the same output file (" + output_path + ").")
                        self._file_dependencies_cumulative[output_path] = copy.deepcopy(input_paths)
                        # now recursively add dependencies of each input file, if any
                        # ATTENTION: Actually, this is an ugly hack, there's no recursion here, but it still
                        # works. It might break at some point, but right now it's fine. TODO: Prove that it
                        # will remain fine.
                        '''
                        p = self.parent
                        l = input_paths
                        for path in l:
                            if path in p._file_dependencies_cumulative:
                                self._file_dependencies_cumulative[output_path].extend(p._file_dependencies_cumulative[path])
                        '''
                            
        # now that the _run_info exists, it remains constant, just return it
        return self._run_info

    def get_run_ids(self):
        return self.get_run_info().keys()

    def get_options_hashtag(self):
        return hashlib.sha1(json.dumps(self.options, sort_keys=True)).hexdigest()[0:4]

    def get_step_id(self):
        return self._step_name

    def get_output_directory(self):
        return os.path.join(self._pipeline.config['destination_path'], 
            '%s-%s' % (self.get_step_id(), self.get_options_hashtag()))

    def get_temp_output_directory(self):
        while True:
            token = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(8))
            path = os.path.join(self._pipeline.config['destination_path'], 'temp', 'temp-' + token)
            if not os.path.exists(path):
                return path

    def get_run_state(self, run_id):

        def path_up_to_date(outpath, inpaths = []):
            if not AbstractStep.fsc.exists(outpath):
                return False
            for inpath in inpaths:
                if not AbstractStep.fsc.exists(inpath):
                    return False
                if AbstractStep.fsc.getmtime(inpath) > AbstractStep.fsc.getmtime(outpath):
                    return False
            return True

        run_info = self.get_run_info()
        all_output_files_exist = True
        all_input_files_exist = True
        for annotation in run_info[run_id]['output_files'].keys():
            for output_file in run_info[run_id]['output_files'][annotation].keys():
                input_files = self._file_dependencies_cumulative[output_file]
                if not path_up_to_date(output_file, input_files):
                    all_output_files_exist = False
                for input_file in input_files:
                    if not path_up_to_date(input_file, []):
                        all_input_files_exist = False

        if all_input_files_exist:
            if all_output_files_exist:
                return self._pipeline.states.FINISHED
            else:
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
        for annotation in temp_run_info['output_files'].keys():
            for out_path, in_paths in temp_run_info['output_files'][annotation].items():
                temp_paths[out_path] = os.path.join(temp_directory, os.path.basename(out_path))

        temp_run_info = fix_dict(temp_run_info, fix_func_dict_subst, temp_paths)

        self._pipeline.notify("[INFO] starting %s/%s" % (str(self), run_id))
        try:
            self.execute(run_id, temp_run_info)
        except Exception as e:
            self._pipeline.notify("[BAD] %s/%s failed: %s" % (str(self), run_id, str(e)))
            raise

        # if we're here, we can assume the step has finished successfully
        # now rename the output files (move from temp directory to
        # destination directory)
        for annotation in temp_run_info['output_files'].keys():
            for out_path in temp_run_info['output_files'][annotation].keys():
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
        
        message = "[OK] %s/%s successfully finished.\n" % (str(self), run_id)
        message += str(self) + ': ' + remaining_task_info + "\n"
        self._pipeline.notify(message)

        # now write the annotation
        log = {}
        log['step'] = {}
        log['step']['options'] = self.options
        log['step']['id'] = self.get_step_id()
        log['run'] = {}
        log['run']['run_info'] = self.get_run_info()[run_id]
        log['run']['run_id'] = run_id
        log['config'] = self._pipeline.config
        log['git_hash_tag'] = self._pipeline.git_hash_tag
        log['tool_versions'] = self._pipeline.tool_versions
        log['pipeline_log'] = unix_pipeline.up_log

        annotation_path = os.path.join(self.get_output_directory(), '.' + run_id + '-annotation.yaml')
        # overwrite the annotation if it already exists
        with open(annotation_path, 'w') as f:
            f.write(yaml.dump(log, default_flow_style = False))

        # create a symbolic link to the annotation for every output file
        for annotation in temp_run_info['output_files'].keys():
            for out_path in temp_run_info['output_files'][annotation].keys():
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

    def tool(self, key):
        '''
        Return full path to a configured tool.
        '''
        return self._pipeline.config['tools'][key]['path']
    
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
