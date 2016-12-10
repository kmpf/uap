#!/usr/bin/env python

import os
import glob
import logging
import abstract_step as abst
import pipeline

logger = logging.getLogger("uap_logger")

class UapFile(object):

    def __init__(self, run, path = None):
        self._os_path_methods = ['abspath', 'basename', 'exists']
        self._path = path
        '''Path to this file (if any)'''
        self._run = run
        '''Run object that created or creats this file'''

    def __getattr__(self, name):
        #nonlocal os_path_methods
        if name in self._os_path_methods:
            return getattr(os.path, name)(self.path)
        else:
            raise StandardError("Could not find method call for %s" % name)
    @property
    def path(self):
        return self._path

    @property
    def run(self):
        return self._run

    def is_source(self):
        '''Returns True if self.run is instance of abst.AbstractSourceStep'''
        return True if self.run.is_source() else False

class UapDataFile(UapFile):

    def __init__(self,
                 run,
                 path = None,
                 upstream_files = list(),
                 downstream_files = list()):
        super(type(self), self).__init__(run, path)
        '''Initialize super object'''
        self._volatile = False
        '''Boolean indicating if this object represents a volatilized file'''

        # Set
        if not self.is_uap_file_list(upstream_files):
            raise StandardError("Expected list of UapDataFile objects")
        self._upstream_files = upstream_files
        '''UapFile objects which preceed this file in the DAG'''
        if not self.is_uap_file_list(downstream_files):
            raise StandardError("Expected list of UapDataFile objects")
        self._downstream_files = downstream_files
        '''UapFile objects which succeed this file in the DAG'''

    @staticmethod
    def is_uap_file_list(uap_file_list):
        '''
        Checks a given list to contain only UapDataFile objects
        '''
        if not type(uap_file_list) is list:
            return False
        for f in uap_file_list:
            if not isinstance(f, UapDataFile):
                return False
        return True
            
    def is_leaf(self):
        '''Returns True if self.downstream_files equals [None]'''
        return True if self.downstream_files == [None] else False
    
    @property
    def downstream_files(self):
        '''
        Returns list of suceeding/downstream files, which depend on
        this file.
        '''
        return self._downstream_files

    @downstream_files.setter
    def downstream_files(self, dfs = list()):
        # Check input dfs
        if not self.is_uap_file_list(dfs):
            raise StandardError("Expected list of UapDataFile objects")
        self._downstream_files = dfs
    
    @property
    def upstream_files(self):
        '''
        Returns list of preceeding/upstream files, which were used to create
        this file.
        '''
        return self._upstream_files

    @upstream_files.setter
    def upstream_files(self, dfs = list()):
        # Check types of upstream files
        if not self.is_uap_file_list(dfs):
            raise StandardError("Expected list of UapDataFile objects")
        self._upstream_files = dfs

    @property
    def volatile(self):
        '''
        Return True, if self.abspath is not found but instead a
        volatile placeholder file exists.
        '''
        return self._volatile
        
class UapPingFile(UapFile):

    def __init__(self, run, path = None):
        super(type(self), self).__init__(run, path)
        

class UapSubmitScriptFile(UapFile):

    def __init__(self, run, path = None):
        super(type(self), self).__init__(run, path)


###############################################################################
# Start doing the essential stuff; fiddling around with file names
###############################################################################

        
def find_files(folder):

    ls = glob.glob(os.path.join(folder, "*"))
    ls.extend(glob.glob(os.path.join(folder, ".*")))
    for i in ls:
        print(os.path.basename(i))
    found_files = [f for f in ls if os.path.isfile(f)]
    found_files += [ i for f in ls for i in find_files(f) if os.path.isdir(f)]
    return [os.path.abspath(f) for f in found_files]

def main(args):
    cwd = os.getcwd()
    '''Current working directory'''

    p = pipeline.Pipeline(arguments=args)
    uap_out_dir = os.path.abspath(p.config["destination_path"])
    '''uap output directory'''

    # Get list of all files managed by current pipeline
    
    ping_files = list()
    submit_script_files = list()
    yaml_files = list()
    data_files = dict()
    '''Dict file path is key; UapDataFile object the value'''
    downstream_files = dict()

    
    # get steps from pipeline
    steps = p.steps
    # get runs from steps (hopefully they know about the files)
    runs = [run for step in steps.values() for run in step.get_runs().values()]
    print("CWD: %s" % cwd)
    # Watch out for the FILES
    for run in runs:
        run_out_dir = run.get_output_directory()
        if run.is_source():
            run_out_dir = ""
            print("Source Step: %s" % run.get_output_files())

        ping_files.append(
            UapPingFile(run, run.get_executing_ping_file()))
        ping_files.append(
            UapPingFile(run, run.get_queued_ping_file()))
        submit_script_files.append(
            UapSubmitScriptFile(run, run.get_submit_script_file()))
        yaml_files
        
        for v in run.get_output_files().values():
            for result_file in v.keys():
                print("** Output file: %s" % result_file)
                out_file = os.path.join(uap_out_dir, run_out_dir, result_file)
                if os.path.exists(out_file):
                    print("** Output file EXISTS!")
                data_files.update({out_file: UapDataFile(run, out_file)})

    # Create the keys in the downstream_files dict
    downstream_files = { abspath: list() for abspath in data_files.keys()}
    for udf in data_files.values():
        if udf.is_source():
            continue
        
        input_files = udf.run.get_input_files_for_output_file(udf.basename)
        # print("%s: %s" % (udf.abspath, input_files))
        # print(udf.run.get_output_directory())

        # Search for upstream UapDataFiles of the current udf (UapDataFile)
        for in_file in input_files:
            upstream_files = list()
            if not os.path.isabs(in_file):
                in_file_abspath = os.path.join(uap_out_dir, in_file)
                # Add UapDataFile to upstream_files list
                upstream_files.append(
                    data_files[in_file_abspath] )
                # Add udf to downstream_files dict for files its directly upstream
                # of
                downstream_files[in_file_abspath].append(udf)
        udf.upstream_files = upstream_files

    # Update downstream_files for all objects in data_files
    for k,v in data_files.items():
        v.downstream_files = downstream_files[k]
    
    all_uap_files = ping_files + submit_script_files + data_files.values()

    # Get list of all files inside uap_out_dir

    all_files = find_files(uap_out_dir)

    # Set operation to exclude files which we know belong to UapFiles
    all_files_set = set(all_files)
    uap_files_set = {i.abspath for i in all_uap_files}
    non_uap_files = all_files_set - uap_files_set
    print('### NOT CAPTURED FILES ###')
    for f in non_uap_files:
        print(os.path.basename(f))
    print('### Ping Files ###')
    for f in ping_files:
        print(f.basename)
