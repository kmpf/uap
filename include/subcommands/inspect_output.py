#!/usr/bin/env python

import os
import glob
import logging

from abstract_step import AbstractSourceStep
import pipeline

logger = logging.getLogger("uap_logger")

class UapFile(object):

    def __init__(self, run, path = None):
        self._os_path_methods = ['abspath', 'exists']
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

        # Check types of upstream and downstream files
        for name, files in {"upstream_files": upstream_files,
                      "downstream_files": downstream_files}.items():
            if not type(files) is list:
                raise StandardError("Parameter %s has to be a list" % name)
            if not files == [None]:
                for f in files:
                    print(type(f))
                    if not isinstance(f, UapDataFile):
                        raise StandardError("Expected UapDataFile got: %s"  % f)
        # Set 
        self._upstream_files = upstream_files
        '''UapFile objects which preceed this file in the DAG'''
        self._downstream_files = downstream_files
        '''UapFile objects which succeed this file in the DAG'''

    @property
    def is_source(self):
        '''Returns True if self.upstream_files equals [None]'''
        return True if self.upstream_files == [None] else False

    @property
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
        
    @property
    def upstream_files(self):
        '''
        Returns list of preceeding/upstream files, which were used to create
        this file.
        '''
        return self._upstream_files

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
    found_files = [f for f in ls if os.path.isfile(f)]
    found_files += [ i for f in ls for i in find_files(f) if os.path.isdir(f)]
    return found_files

def main(args):
    cwd = os.getcwd()
    '''Current working directory'''

    p = pipeline.Pipeline(arguments=args)
    output_dir = os.path.abspath(p.config["destination_path"])
    '''uap output directory'''

    # Get list of all files managed by current pipeline
    
    ping_files = list()
    submit_script_files = list()
    output_files = list()
    
    # get steps from pipeline
    steps = p.steps
    # get runs from steps (hopefully they know about the files)
    runs = [run for step in steps.values() for run in step.get_runs().values()]

#     for k,v in p.input_files_for_task_id.items():
#         for f in v:
#             print("%s: %s" % (k, f))


        
    # Watch out for the FILES
    for run in runs:
        if isinstance(run, AbstractSourceStep):
            print(run.get_output_files())
        run_out_dir = run.get_output_directory()
        ping_files.append(UapPingFile(run, run.get_executing_ping_file()))
        ping_files.append(UapPingFile(run, run.get_queued_ping_file()))
        submit_script_files.append(
            UapSubmitScriptFile(
                run,
                run.get_submit_script_file()
            )
        )
        for v in run.get_output_files().values():
            print("== %s ==" % run.get_run_id())
            for result_file in v.keys():
                out_file = os.path.join(cwd, run_out_dir, result_file)
                s = out_file
                if os.path.exists(out_file):
                    s += " -> Missing!"
                print(s)
                output_files.append(UapDataFile(run, out_file))


#        print(output_files)

    all_uap_files = ping_files + submit_script_files + output_files
#    all_uap_files = [ p.abspath() for p in all_uap_files]
    for i in all_uap_files:
        #print(i.path)
        print(i.abspath)
        print(i.exists)

    uap_files_set = set(all_uap_files)
        
#    print("all_uap_files %s" % len(all_uap_files))

    # Get list of all files inside output_dir

    all_files = find_files(output_dir)
#     for i in all_files:
#         print(i)

    exit(1)
    all_files_set = set(all_files)
    non_uap_files = all_files_set - uap_files_set
    for f in non_uap_files:
        print(f)
#    print("all_files length %s" % len(all_files))
