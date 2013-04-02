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
import traceback
import yaml

def get_step_class_for_key(key):
    classes = [_ for _ in inspect.getmembers(__import__(key), inspect.isclass) if AbstractStep in _[1].__bases__]
    if len(classes) != 1:
        raise StandardError("need exactly one subclass of AbstractStep in " + key)
    return classes[0][1]

'''
Turns filenames like these (paths):
    ['RIB0000794-cutadapt-R1.fastq.gz', 'RIB0000794-cutadapt-R2.fastq.gz']
and this (tags):
    ['R1', 'R2']
...into this:
    {
        'R1': 'RIB0000794-cutadapt-R1.fastq.gz',
        'R2': 'RIB0000794-cutadapt-R2.fastq.gz'
    }
If this is not possible without ambiguities, a StandardError is thrown
'''
def assign_strings(paths, tags):

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

class AbstractStep(object):
    def __init__(self, pipeline):
        self.dependencies = []
        self.options = {}
        self.pipeline = pipeline
        self.step_name = self.__module__
        self._run_info = None
        self._cores = 1

    def set_cores(self, cores):
        self._cores = cores

    def set_options(self, options):
        self.options = options

    def add_dependency(self, parent):
        if not isinstance(parent, AbstractStep):
            raise StandardError("parent argument must be an AbstractStep")
        if parent == self:
            raise StandardError("cannot add a node as its own dependency")
        self.dependencies.append(parent)

    def get_input_run_info(self):
        if len(self.dependencies) == 0:
            raise StandardError("You asked for input files of a step with no dependencies. This shouldn't happen.")
        elif len(self.dependencies) == 1:
            return copy.deepcopy(self.dependencies[0].get_run_info())
        else:
            raise NotImplementedError("DAG not implemented yet.")

    def setup_runs(self, input_run_info):
        raise NotImplementedError()

    def get_run_info(self):
        # create run info if it doesn't exist yet
        if self._run_info == None:
            input_run_info = None
            full_paths = dict()
            if len(self.dependencies) > 0:
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

            if len(self.dependencies) > 0:

                for run_id, _ in self._run_info.items():
                    for annotation in self._run_info[run_id]['output_files'].keys():
                        for path in self._run_info[run_id]['output_files'][annotation].keys():
                            full_paths[path] = os.path.join(self.get_output_directory(), path)

                self._run_info = fix_dict(self._run_info, fix_func_dict_subst, full_paths)

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
        while len(p.dependencies) > 0:
            if len(p.dependencies) > 1:
                raise NotImplementedError("Full DAG not implemented yet, trees only for now.")
            p = p.dependencies[0]
            if p.__module__ != 'source':
                path.append(p.__module__)
                if with_options:
                    path[-1] += '-' + p.get_options_hashtag()
        path.reverse()
        return path

    def get_options_hashtag(self):
        return hashlib.sha1(json.dumps(self.options, sort_keys=True)).hexdigest()[0:8]

    def get_step_id(self):
        return '/'.join(self.get_dependency_path())

    def get_output_directory(self):
        # add the options to the output path so that different options result in
        # a different directory
        dependency_path = self.get_dependency_path(True)
        if (self.pipeline.run_mode == self.pipeline.run_modes.TEST_RUN):
            dependency_path.insert(0, 'test')
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
                if not os.path.exists(outpath):
                    return False
                for inpath in inpaths:
                    if not os.path.exists(inpath):
                        return False
                    if os.path.getmtime(inpath) > os.path.getmtime(outpath):
                        return False
                return True

        run_info = self.get_run_info()
        all_output_files_exist = True
        all_input_files_exist = True
        for annotation in run_info[run_id]['output_files'].keys():
            for output_file, input_files in run_info[run_id]['output_files'][annotation].items():
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
            os.makedirs(temp_directory)

            # call execute() but pass output file paths with the temporary directory
            temp_run_info = copy.deepcopy(self.get_run_info()[run_id])
            temp_paths = {}
            for annotation in temp_run_info['output_files'].keys():
                for out_path, in_paths in temp_run_info['output_files'][annotation].items():
                    temp_paths[out_path] = os.path.join(temp_directory, os.path.basename(out_path))

            temp_run_info = fix_dict(temp_run_info, fix_func_dict_subst, temp_paths)

            self.pipeline.notify("[INFO] starting " + self.get_step_id() + "/" + run_id)
            try:
                self.execute(run_id, temp_run_info)
            except Exception as e:
                self.pipeline.notify("[BAD] " + self.get_step_id() + "/" + run_id + " failed: " + str(e))
                raise

            self.pipeline.notify("[OK] " + self.get_step_id() + "/" + run_id + " successfully finished")

            # if we're here, we can assume the step has finished successfully
            # now rename the output files (move from temp directory to
            # destination directory)
            for annotation in temp_run_info['output_files'].keys():
                for out_path in temp_run_info['output_files'][annotation].keys():
                    destination_path = os.path.join(self.get_output_directory(), os.path.basename(out_path))
                    # TODO: if the destination path already exists, this will overwrite the file.
                    os.rename(out_path, destination_path)

            # now write the annotation
            log = {}
            log['all_samples'] = self.pipeline.all_samples
            log['step'] = {}
            log['step']['options'] = self.options
            log['step']['id'] = self.get_step_id()
            log['run_info'] = self.get_run_info()
            log['config'] = self.pipeline.config
            log['git_hash_tag'] = self.pipeline.git_hash_tag

            annotation_path = os.path.join(self.get_output_directory(), 'annotation-' + hashlib.sha1(json.dumps(log, sort_keys=True)).hexdigest()[0:8] + '.yaml')
            with open(annotation_path, 'w') as f:
                f.write(yaml.dump(log, default_flow_style = False))

            # finally, remove the temporary directory if it's empty
            try:
                os.rmdir(temp_directory)
            except OSError:
                pass

    def execute(self, run_id, run_info):
        raise NotImplementedError()

    def tool(self, key):
        return self.pipeline.config['tools'][key]['path']

    def __str__(self):
        s = self.step_name
        #if len(self.options.keys()) > 0:
            #s += " with options: " + str(self.options)
        #if len(self.dependencies) > 0:
            #s += " with dependencies: " + ', '.join(['(' + str(_) + ')' for _ in self.dependencies])
        return s

