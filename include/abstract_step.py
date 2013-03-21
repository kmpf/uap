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
import yaml

def get_step_class_for_key(key):
    classes = [_ for _ in inspect.getmembers(__import__(key), inspect.isclass) if AbstractStep in _[1].__bases__]
    if len(classes) != 1:
        raise StandardError("need exactly one subclass of AbstractStep in " + key)
    return classes[0][1]

class AbstractStep(object):
    def __init__(self, pipeline):
        self.dependencies = []
        self.options = {}
        self.pipeline = pipeline
        self.step_name = self.__module__
        self.run_info = None

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
            return copy.copy(self.dependencies[0].get_run_info())
        else:
            raise NotImplementedError("DAG not implemented yet.")

    def setup_runs(self, input_run_info):
        raise NotImplementedError()

    def get_run_info(self):
        # create run info if it doesn't exist yet
        if self.run_info == None:
            input_run_info = None
            full_paths = dict()
            if len(self.dependencies) > 0:
                # it's not the source step
                input_run_info = self.get_input_run_info()
                # strip directories from file names, strip input files
                for run_id, output_files in input_run_info.items():
                    for path in output_files.keys():
                        basename = os.path.basename(path)
                        if basename in full_paths:
                            raise StandardError("There are multiple input filenames with the same basename.")
                        full_paths[basename] = path
                    input_run_info[run_id] = sorted([os.path.basename(path) for path in output_files.keys()])

            self.run_info = self.setup_runs(input_run_info)
            if len(self.dependencies) > 0:
                run_info_temp = self.run_info
                self.run_info = {}
                # re-add directories to input and output file names
                for run_id, output_files in run_info_temp.items():
                    self.run_info[run_id] = {}
                    for path in output_files.keys():
                        full_path = os.path.join(self.get_output_directory(), path)
                        self.run_info[run_id][full_path] = []
                        for in_path in output_files[path]:
                            self.run_info[run_id][full_path].append(full_paths[in_path])
        return self.run_info

    def get_run_ids(self):
        run_info = self.get_run_info()
        return run_info.keys()

    def get_dependency_path(self):
        path = []
        p = self
        path.append(p.__module__)
        while len(p.dependencies) > 0:
            if len(p.dependencies) > 1:
                raise NotImplementedError("Full DAG not implemented yet, trees only for now.")
            p = p.dependencies[0]
            if p.__module__ != 'source':
                path.append(p.__module__)
        path.reverse()
        return path

    def get_step_id(self):
        return '/'.join(self.get_dependency_path())

    def get_output_directory(self):
        # add the options to the output path so that different options result in
        # a different directory
        options_checksum = hashlib.sha1(json.dumps(self.options, sort_keys=True)).hexdigest()[0:8]
        dependency_path = self.get_dependency_path()
        if (self.pipeline.run_mode == pipeline.Pipeline.TEST_RUN):
            dependency_path[-1] += '-test-run'
        dependency_path[-1] += '-' + options_checksum
        return os.path.join(self.pipeline.config['destinationPath'], *dependency_path)

    def get_temp_output_directory(self):
        while True:
            token = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(8))
            path = os.path.join(self.pipeline.config['destinationPath'], 'temp', 'temp-' + token)
            if not os.path.exists(path):
                return path

    def get_run_state(self, run_id, dry_run_cache = None):

        def path_up_to_date(outpath, inpaths = [], dry_run_cache = None):
            if dry_run_cache == None:
                # check the file system
                if not os.path.exists(outpath):
                    return False
                for inpath in inpaths:
                    if not os.path.exists(inpath):
                        return False
                    if os.path.getmtime(inpath) > os.path.getmtime(outpath):
                        return False
                return True
            else:
                # this is a dry run, use the dry run cache
                if not outpath in dry_run_cache:
                    return False
                for inpath in inpaths:
                    if not inpath in dry_run_cache:
                        return False
                    if dry_run_cache[inpath] > dry_run_cache[outpath]:
                        return False
                return True

        run_info = self.get_run_info()
        all_output_files_exist = True
        all_input_files_exist = True
        # TODO: check whether output files are up-to-date regarding their input files
        for output_file, input_files in run_info[run_id].items():
            if not path_up_to_date(output_file, input_files, dry_run_cache):
                all_output_files_exist = False
            for input_file in input_files:
                if not path_up_to_date(input_file, [], dry_run_cache):
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
        # create the output directory if it doesn't exist yet
        if not os.path.isdir(self.get_output_directory()):
            os.makedirs(self.get_output_directory())
        # also create a temporary directory for the output file
        temp_directory = self.get_temp_output_directory()
        os.makedirs(temp_directory)

        # call execute() but pass output file paths with the temporary directory
        temp_run_info = {}
        for out_path, in_paths in self.get_run_info()[run_id].items():
            temp_out_path = os.path.join(temp_directory, os.path.basename(out_path))
            temp_run_info[temp_out_path] = in_paths

        print("executing " + self.get_step_id() + " " + run_id)
        self.execute(run_id, temp_run_info)

        # if we're here, we can assume the step has finished successfully
        # now rename the output files (move from temp directory to
        # destination directory)
        for out_path in temp_run_info.keys():
            destination_path = os.path.join(self.get_output_directory(), os.path.basename(out_path))
            os.rename(out_path, destination_path)

        # finally, remove the temporary directory if it's empty
        try:
            os.rmdir(temp_directory)
        except OSError:
            pass

    def execute(self, run_id, run_info):
        print("WARNING: Just creating empty output files for " + self.get_step_id() + "/" + run_id + " due to missing implementation.")
        for path in run_info.keys():
            with open(path, 'w') as f:
                pass

    def __str__(self):
        s = self.step_name
        #if len(self.options.keys()) > 0:
            #s += " with options: " + str(self.options)
        #if len(self.dependencies) > 0:
            #s += " with dependencies: " + ', '.join(['(' + str(_) + ')' for _ in self.dependencies])
        return s

