import sys
sys.path.append('./include/steps')
import copy
import datetime
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

def get_step_class_for_key(key):
    classes = [_ for _ in inspect.getmembers(__import__(key), inspect.isclass) if AbstractStep in _[1].__bases__]
    if len(classes) != 1:
        raise StandardError("need exactly one subclass of AbstractStep in " + key)
    return classes[0][1]

def assign_strings(paths, tags):
    '''
    Assign strings (path names, for example) to tags. Example:
    
      - paths = ['RIB0000794-cutadapt-R1.fastq.gz', 'RIB0000794-cutadapt-R2.fastq.gz']
      - tags = ['R1', 'R2']
      - result = { 'R1': 'RIB0000794-cutadapt-R1.fastq.gz', 'R2': 'RIB0000794-cutadapt-R2.fastq.gz' }
      
    If this is not possible without ambiguities, a StandardError is thrown.
    '''

    def check_candidate(paths, tags, head, tail):
        chopped = []
        for path in paths:
            if path[:len(head)] != head:
                return None
            if path[-len(tail):] != tail:
                return None
            chopped.append((path[len(head):-len(tail)], path))

        if [_[0] for _ in sorted(chopped)] == sorted(tags):
            result = {}
            for _ in sorted(chopped):
                result[_[0]] = _[1]
            return result

        return None


    results = {}
    if len(paths) != len(tags):
        raise StandardError("Number of tags must be equal to number of paths")
    for tag in tags:
        for path in paths:
            result_candidate = {}
            if tag in path:
                # find all occurences of tag in path
                offset = 0
                while path.find(tag, offset) >= 0:
                    index = path.find(tag, offset)
                    head = path[:index]
                    tail = path[(index+len(tag)):]
                    # now try chopping off head and tail from every path
                    # and see whether we can unambiguously assign a path
                    # to every tag, if yes, we have a result candidate
                    result_candidate = check_candidate(paths, tags, head, tail)
                    if result_candidate:
                        results[json.dumps(result_candidate, sort_keys = True)] = result_candidate
                    offset = index + 1
    if len(results) != 1:
        raise StandardError("Unable to find an unambiguous mapping.")
    return results[results.keys()[0]]

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

# fs_cache records return values of 'exists' and 'getmtime' which is done
# so that actual syscalls for this are only executed once and later re-used
# from the cache

fs_cache = dict()

def fs_cache_path_exists(path):
    key = 'exists/' + path
    if key in fs_cache:
        return fs_cache[key]
    result = os.path.exists(path)
    fs_cache[key] = result
    return result

def fs_cache_path_getmtime(path):
    key = 'getmtime/' + path
    if key in fs_cache:
        return fs_cache[key]
    result = os.path.getmtime(path)
    fs_cache[key] = result
    return result

def fs_cache_flush():
    fs_cache.clear()

class AbstractStep(object):
    def __init__(self, pipeline):
        self.parent = None
        self.options = {}
        self.pipeline = pipeline
        self.step_name = self.__module__
        self._run_info = None
        # _file_dependencies_cumulative is secondary data which gets generated 
        # from _run_info and keeps track of each output file's file 
        # dependencies (regardless of the run id) and by file dependencies
        # I mean ALL file dependencies including those from all of its
        # parent steps (yup, that may become a lot of files)
        self._file_dependencies_cumulative = {}
        self._cores = 1
        self._temp_directory = None

    def set_cores(self, cores):
        self._cores = cores

    def set_options(self, options):
        self.options = options

    def add_dependency(self, parent):
        if not isinstance(parent, AbstractStep):
            raise StandardError("Error: parent argument must be an AbstractStep.")
        if parent == self:
            raise StandardError("Cannot add a node as its own dependency.")
        if self.parent:
            raise StandardError("This step already has a parent.")
        self.parent = parent

    def get_input_run_info(self):
        if not self.parent:
            raise StandardError("You asked for the input run info of a step with no parent. This shouldn't happen.")
        return copy.deepcopy(self.parent.get_run_info())

    def setup_runs(self, input_run_info):
        raise NotImplementedError()

    def get_run_info(self):
        # create run info if it doesn't exist yet
        if self._run_info == None:
            # create input run info and simplify it a bit for setup_runs
            input_run_info = None
            full_paths = dict()
            if self.parent:
                # it's not the source step
                input_run_info = copy.deepcopy(self.get_input_run_info())

                # strip directories from file names, strip input files
                for run_id, run_info in input_run_info.items():
                    for annotation in run_info['output_files'].keys():
                        for path in run_info['output_files'][annotation].keys():
                            basename = os.path.basename(path)
                            if basename in full_paths:
                                raise StandardError("There are multiple input filenames with the same basename.")
                            full_paths[basename] = path
                            run_info['output_files'][annotation][basename] = run_info['output_files'][annotation][path]
                            run_info['output_files'][annotation].pop(path)

            self._run_info = self.setup_runs(input_run_info)
            for run_id in self._run_info.keys():
                if '/' in run_id:
                    raise StandardError("A run_id must not contain a '/': " + run_id + " for " + str(self))

            if self.parent:

                for run_id, _ in self._run_info.items():
                    for annotation in self._run_info[run_id]['output_files'].keys():
                        for path in self._run_info[run_id]['output_files'][annotation].keys():
                            full_paths[path] = os.path.join(self.get_output_directory(), path)

                self._run_info = fix_dict(self._run_info, fix_func_dict_subst, full_paths)
                
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
                            p = self.parent
                            l = input_paths
                            for path in l:
                                if path in p._file_dependencies_cumulative:
                                    self._file_dependencies_cumulative[output_path].extend(p._file_dependencies_cumulative[path])
                            
        # now that the _run_info exists, it remains constant, just return it
        return self._run_info

    def get_run_ids(self):
        run_info = self.get_run_info()
        return run_info.keys()

    def get_dependency_path(self, with_options = False):
        path = []
        p = self
        path.append(p.__module__)
        if with_options:
            path[-1] += '-' + p.get_options_hashtag()
        while p.parent:
            p = p.parent
            if p.__module__ != 'source':
                path.append(p.__module__)
                if with_options:
                    path[-1] += '-' + p.get_options_hashtag()
        path.reverse()
        return path

    def get_options_hashtag(self):
        return hashlib.sha1(json.dumps(self.options, sort_keys=True)).hexdigest()[0:4]

    def get_step_id(self):
        return '/'.join(self.get_dependency_path())

    def get_task_id_fragmented(self, run_id):
        '''
        get_task_id_fragmented returns the step ID, but in a way which puts
        the most important information first, e. g.
        ``segemehl/Sample_COPD_363 (cutadapt/fix_cutadapt/...)``
        instead of
        ``cutadapt/fix_cutadapt/segemehl/Sample_COPD_363``
        '''
        path = self.get_dependency_path()
        result = path[-1] + '/' + run_id
        if len(path) > 1:
            result += ' (' + '/'.join(path[:-1]) + '/...)'
        return result

    def get_output_directory(self):
        # add the options to the output path so that different options result in
        # a different directory
        dependency_path = self.get_dependency_path(True)
        return os.path.join(self.pipeline.config['destination_path'], *dependency_path)

    def get_temp_output_directory(self):
        while True:
            token = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(8))
            path = os.path.join(self.pipeline.config['destination_path'], 'temp', 'temp-' + token)
            if not os.path.exists(path):
                return path

    def get_run_state(self, run_id):

        def path_up_to_date(outpath, inpaths = []):
            if self.pipeline.run_mode == self.pipeline.run_modes.DRY_RUN:
                # this is a dry run, use the dry run cache
                if not outpath in self.pipeline.dry_run_cache:
                    return False
                for inpath in inpaths:
                    if not inpath in self.pipeline.dry_run_cache:
                        return False
                    if self.pipeline.dry_run_cache[inpath] > self.pipeline.dry_run_cache[outpath]:
                        return False
                return True
            else:
                # check the file system
                if not fs_cache_path_exists(outpath):
                    return False
                for inpath in inpaths:
                    if not fs_cache_path_exists(inpath):
                        return False
                    if fs_cache_path_getmtime(inpath) > fs_cache_path_getmtime(outpath):
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
                return self.pipeline.states.FINISHED
            else:
                return self.pipeline.states.READY
        else:
            return self.pipeline.states.WAITING

    def dry_run(self, run_id, dry_run_cache):
        run_info = self.get_run_info()[run_id]
        for path in run_info.keys():
            dry_run_cache[path] = datetime.datetime.now()

    def run(self, run_id):
        if self.pipeline.run_mode == self.pipeline.run_modes.DRY_RUN:
            print("dry-running " + self.get_step_id() + "/" + run_id)
            run_info = copy.deepcopy(self.get_run_info()[run_id])
            for annotation in run_info['output_files']:
                for out_path in run_info['output_files'][annotation].keys():
                    self.pipeline.dry_run_cache[out_path] = datetime.datetime.now()
        else:
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

            self.pipeline.notify("[INFO] starting " + self.get_task_id_fragmented(run_id))
            try:
                self.execute(run_id, temp_run_info)
            except Exception as e:
                self.pipeline.notify("[BAD] " + self.get_task_id_fragmented(run_id) + " failed: " + str(e))
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
            fs_cache_flush()
            
            count = {}
            for _ in self.get_run_ids():
                state = self.get_run_state(_)
                if not state in count:
                    count[state] = 0
                count[state] += 1
            remaining_task_info = ', '.join([str(count[_]) + ' ' + _.lower() for _ in sorted(count.keys())])
            
            message = "[OK] " + self.get_task_id_fragmented(run_id) + " successfully finished.\n"
            message += str(self) + ': ' + remaining_task_info + "\n"
            self.pipeline.notify(message)

            # now write the annotation
            log = {}
            log['all_samples'] = self.pipeline.all_samples
            log['step'] = {}
            log['step']['options'] = self.options
            log['step']['id'] = self.get_step_id()
            log['run'] = {}
            log['run']['run_info'] = self.get_run_info()[run_id]
            log['run']['run_id'] = run_id
            log['config'] = self.pipeline.config
            log['git_hash_tag'] = self.pipeline.git_hash_tag
            log['tool_versions'] = self.pipeline.tool_versions
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

    def execute(self, run_id, run_info):
        raise NotImplementedError()

    def tool(self, key):
        return self.pipeline.config['tools'][key]['path']
    
    def get_temporary_path(self, prefix, suffix):
        '''
        Returns a temporary path with the prefix and suffix specified. The
        returned path will be in the temporary directory.
        '''
        if not self._temp_directory:
            raise StandardError("Temporary directory not set, you cannot call get_temporary_path from setup_runs.")

        _, _path = tempfile.mkstemp(suffix, prefix, self._temp_directory)
        os.close(_)
        os.unlink(_path)

        return _path        

    def __str__(self):
        return self.step_name

