import sys
from datetime import datetime
import glob
import json
from logging import getLogger
import os
import random
import stat
import string
import tempfile

import yaml

import abstract_step as abst
import command as command_info
import exec_group
import pipeline_info
import misc

logger = getLogger("uap_logger")

class Run(object):
    '''
    The Run class is a helper class which represents a run in a step. Declare
    runs inside AbstractStep.runs() via::
    
        with self.new_run(run_id) as run:
            # declare output files, private and public info here
            
    After that, use the available methods to configure the run.
    The run has typically no information about input connections only about
    input files.
    '''
    def __init__(self, step, run_id):
        if '/' in run_id:
            logger.error("Error: A run ID must not contain a slash: %s." % run_id)
            sys.exit(1)
        self._step = step
        '''
        Step this run belongs to.
        '''
        self._run_id = run_id
        '''
        Identifier of this run.
        '''
        self._private_info = dict()
        self._public_info = dict()
        self._input_files = list()
        self._output_files = dict()
        for out_connection in self._step.get_out_connections():
            self.add_out_connection(out_connection)
        '''
        Dictionary containing the output files for each outgoing connection and
        their corresponding input files::

           annotation_1:
               out_path_1: [in_path_1, in_path_2, ...]
               out_path_2: ...
           annotation_2: ...

        '''
        self._ping_files = {
            'run': None,
            'queued': None
        }
        self._submit_script = None
        self._exec_groups = list()
        self._temp_paths = list()
        '''
        List of temporary paths which can be either files or paths
        '''
        self._temp_directory = None
        '''
        Contains path to currently used temporary directory if set.
        '''
        self._known_paths = dict()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def new_exec_group(self):
        eg = exec_group.ExecGroup(self)
        self._exec_groups.append(eg)
        return eg
        
    def get_exec_groups(self):
        return self._exec_groups

    def get_step(self):
        return self._step

    def set_run_id(self, run_id):
        self._run_id = run_id
        
    def get_run_id(self):
        return self._run_id

    def get_out_connections(self):
        return self._output_files.keys()

    def get_out_connection(self, connection):
        if not connection.startswith("out/"):
            connection = 'out/' + connection
        if connection in self.get_out_connections():
            return connection
        else:
            logger.error("Connection %s not declared for step %s" %
                         (connection, self.get_step()))
            sys.exit(1)

    def _get_ping_file(self, key):
        if self._ping_files[key] == None:
            self._ping_files[key] = os.path.join(
                self.get_output_directory(),
                '.%s-%s-ping.yaml' % (self.get_run_id(), key)
            )
        return self._ping_files[key]

    def get_executing_ping_file(self):
        return self._get_ping_file('run')

    def get_queued_ping_file(self):
        return self._get_ping_file('queued')

    def get_submit_script_file(self):
        if self._submit_scripts == None:
            self._submit_scripts = os.path.join(
                self.get_output_directory(),
                ".submit-%s-%s.sh" % (self.get_step().get_step_name(),
                                      self.get_run_id())
            )
        return self._submit_scripts

    def replace_output_dir_du_jour(func):
        def inner(self, *args, **kwargs):
            # Collect info to replace du_jour placeholder with temp_out_dir
            step = self.get_step()
            placeholder = self.get_output_directory_du_jour_placeholder()
            temp_out_dir = self.get_output_directory_du_jour()
            
            value = None
            ret_value = func(self, *args, **kwargs)
            # If currently calling AbstractStep.runs() do nothing
            if temp_out_dir == None:
                value = ret_value
            elif isinstance(ret_value, list):
                value = list()
                for string in ret_value:
                    if string != None and placeholder in string:
                        value.append(
                            string.replace(placeholder, temp_out_dir))
                    else:
                        value.append(string)
            elif isinstance(ret_value, str):
                if ret_value != None and placeholder in ret_value:
                    value = ret_value.replace(placeholder, temp_out_dir)
            elif isinstance(ret_value, dict):
                for key in ret_value.keys():
                    if key != None and placeholder in key:
                        new_key = key.replace(placeholder, temp_out_dir)
                        ret_value[new_key] = ret_value.pop(key)
                value = ret_value
            elif ret_value == None:
                value = None
            else:
                logger.error("Function %s does not return list or string object"
                             % func.__class__.__name__)
                sys.exit(1)
            return(value)
        return(inner)

    @replace_output_dir_du_jour
    def get_known_paths(self):
        return self._known_paths

    def add_known_paths(self, known_paths_dict):
        self._known_paths.update(known_paths_dict)

    @replace_output_dir_du_jour
    def get_temp_paths(self):
        '''
        Returns a list of all temporary paths which belong to this run.
        '''
        return self._temp_paths

    def get_output_directory_du_jour_placeholder(self):
        '''
        Returns a placeholder for the temporary output directory, which
        needs to be replaced by the actual temp directory inside the
        abstract_step.execute() method
        '''
        return("<%s-output-directory-du-jour>" %  
               str(self.get_step().__class__.__name__))

    def get_output_directory_du_jour(self):
        '''
        Returns the state-dependent output directory of the step.


        Returns this steps output directory according to its current
        state:
            - if we are currently calling a step's declare_runs()
              method, this will return None
            - if we are currently calling a step's execute() method,
              this will return the temporary directory
            - otherwise, it will return the real output directory
        '''
        if self.get_step()._state == abst.AbstractStep.states.DEFAULT:
            return self.get_output_directory()
        elif self.get_step()._state == abst.AbstractStep.states.EXECUTING:
            return self.get_temp_output_directory()
        else:
            return None

    def get_temp_output_directory(self):
        '''
        Returns the temporary output directory of a run.
        '''
        if self._temp_directory == None:
            while True:
                current_time = datetime.now().strftime('%y%m%d-%H%M%S-%f')
                path = os.path.join(
                    self.get_step().get_pipeline().config['destination_path'],
                    'temp',
                    'temp-%s-%s-%s' % (self.get_step().get_step_name(),
                                       self.get_run_id(), current_time))
                if not os.path.exists(path):
                    self._temp_directory = path
                    return self._temp_directory
        
        return self._temp_directory

    def get_execution_hashtag(self):
        '''
        Creates a hash tag based on the commands to be executed.

        This causes runs to be marked for rerunning if the commands to be
        executed change.
        '''

        # Store step state
        previous_state = self.get_step()._state
        # Set step state to DECLARING to avoid circular dependencies
        self.get_step()._state = abst.AbstractStep.states.DECLARING

        cmd_by_eg = dict()
        eg_count = 0
        for exec_group in self.get_exec_groups():
            eg_count += 1
            cmd_by_eg[eg_count] = dict()
            pipe_count, cmd_count = (0, 0)
            for poc in exec_group.get_pipes_and_commands():
                # for each pipe or command (poc)
                # check if it is a pipeline ...
                if isinstance(poc, pipeline_info.PipelineInfo):
                    pipe_count += 1
                    cmd_by_eg[eg_count]['Pipe %s' % pipe_count] = list()
                    for command in poc.get_commands():
                        cmd_by_eg[eg_count]['Pipe %s' % pipe_count].append(
                            command.get_command())
                # ... or a command
                elif isinstance(poc, command_info.CommandInfo):
                    cmd_count += 1
                    cmd_by_eg[eg_count]['Cmd %s' % cmd_count] = poc.get_command()

        # Set step state back to original state
        self.get_step()._state = previous_state
        return misc.str_to_sha1_b62(json.dumps(cmd_by_eg))[0:8]

    def get_output_directory(self):
        '''
        Returns the final output directory.
        '''
        return os.path.join(
            self.get_step().get_pipeline().config['destination_path'], 
            self.get_step().get_step_name(),
            '%s-%s' % (self.get_run_id(), self.get_execution_hashtag())
        )

    def get_basic_state(self):
        '''
        Determines basic run state of a run.

        Determine the basic run state of a run, which is, at any time, one of
        **waiting**, **ready**, or **finished**.
        
        These states are determined from the current configuration and the
        timestamps of result files present in the file system. In addition to
        these three basic states, there are two additional states which are
        less reliable (see *get_run_state()*).
        '''

        def volatile_path_good(volatile_path, recurse = True):
            '''
            This function receives a volatile path and tries to load the
            placeholder YAML data structure. It then checks all downstream
            paths, which may in turn be volatile placeholder files.
            '''
            
            # reconstruct original path from volatile placeholder path
            path = volatile_path[:-len(AbstractStep.VOLATILE_SUFFIX)]
            
            if AbstractStep.fsc.exists(path):
                # the original file still exists, ignore volatile placeholder
                return False
            
            if not path in self.get_step().get_pipeline()\
                                          .task_id_for_output_file:
                # there is no task which creates the output file
                return False
            
            task_id = self.get_step().get_pipeline().task_id_for_output_file[path]
            
            task = self.get_step().get_pipeline().task_for_task_id[task_id]
#            if not task.step.options['_volatile']:
            if not task.step._options['_volatile']:
                # the task is not declared volatile
                return False
            
            if not AbstractStep.fsc.exists(volatile_path):
                # the volatile placeholder does not exist
                return False
            
            if not recurse:
                return True
            
            try:
                # try to parse the YAML contents
                info = AbstractStep.fsc.load_yaml_from_file(volatile_path)
            except yaml.scanner.ScannerError:
                # error scanning YAML
                return False
            
            # now check whether all downstream files are in place and up-to-date
            # also check whether all downstream files as defined in
            # file_dependencies_reverse are covered

            uncovered_files = set()
            if path in self.get_step().get_pipeline().file_dependencies_reverse:
                uncovered_files = self.get_step().get_pipeline()\
                                                 .file_dependencies_reverse[path]
                
            for downstream_path, downstream_info in info['downstream'].items():
                if downstream_path in self.get_step().get_pipeline()\
                                                     .task_id_for_output_file:
                    # only check this downstream file if there's a task which 
                    # creates it, otherwise, it may be a file which is no more
                    # used
                    pv_downstream_path = change_to_volatile_if_need_be(
                        downstream_path, recurse = False)
                    if not AbstractStep.fsc.exists(pv_downstream_path):
                        return False
                    if not AbstractStep.fsc.getmtime(pv_downstream_path) >= \
                       info['self']['mtime']:
                        return False
                    if downstream_path in uncovered_files:
                        uncovered_files.remove(downstream_path)
                
            if len(uncovered_files) > 0:
                # there are still files defined which are not covered by the
                # placeholder
                return False
                
            return True
        
        def change_to_volatile_if_need_be(path, recurse = True):
            '''
            Changes the file path to volatile path if necessary.
            '''
            if path != None:
                if not AbstractStep.fsc.exists(path):
                    # the real output file does not exist
                    volatile_path = path + AbstractStep.VOLATILE_SUFFIX
                    if volatile_path_good(volatile_path, recurse):
                        return volatile_path
                return path

        def is_path_up_to_date(outpath, inpaths):
            '''
            First, replace paths with volatile paths if the step is marked
            as volatile and the real path is missing.
            But, only consider volatile placeholders if all child tasks
            are finished. That means if a child of a volatile
            step needs to be run because it has been added or an existing step
            has been modified, the volatile placeholders are ignored, thus
            turning the task from 'finished' to 'ready' or 'waiting'
            Hint: The pv_ prefix is for 'possibly volatile'
            '''
            pv_outpath = outpath
            pv_inpaths = list()
            
            if outpath in self.get_step().get_pipeline().task_id_for_output_file:
                pv_outpath = change_to_volatile_if_need_be(outpath)
                
            for inpath in inpaths:
                pv_inpaths.append(change_to_volatile_if_need_be(inpath))
                
            if not AbstractStep.fsc.exists(pv_outpath):
                return False
            for pv_inpath in pv_inpaths:
                if not AbstractStep.fsc.exists(pv_inpath):
                    return False
                if AbstractStep.fsc.getmtime(pv_inpath) > \
                   AbstractStep.fsc.getmtime(pv_outpath):
                    return False
            return True
            
        def up_to_dateness_level(path, level = 0):
            result = level
            if path != None:
                dep_paths = self.get_step().get_pipeline().file_dependencies[path]
                if not is_path_up_to_date(path, dep_paths):
                    result = level + 1
                for dep_path in dep_paths:
                    recursive_result = up_to_dateness_level(dep_path, level + 1)
                    if recursive_result > level + 1:
                        result = max(result, recursive_result)
                return result

        '''
        - finished: all output files exist AND up to date (recursively)
        - ready: NOT all output files exist AND all input files exist AND up to
                 date (recursively)
        - waiting: otherwise
        - if it's ready, it might be executing or queued -> check execute and
          queue ping
        - if it's waiting, it might be queued -> check queue ping
        
        the ping works like this (this example is for execute, same goes for 
        queued):
          - there's a ping file for every task ( = step + run)
          - it contains information about when, how, where the job was started
            etc.
          - its timestamp gets renewed every 30 seconds (touch)
          - as soon as the job has finished, the execute ping file is removed,
            this should also work if the job crashes (however, it cannot work
            if the controlling script receives SIGKILL
          - if its timestamp is no more than 5 minutes old, it is regarded as
            currently executing
          - otherwise, a warning is printed because the ping file is probably 
            stale (no automatic cleanup is performed, manual intervention is
            necessary)
          - warning: this requires all involved systems or the file system to
            be time-synchronized
        '''
        
        run_info = self.get_runs()
        max_level = 0
        for tag, output_files in self.get_output_files_abspath()\
                                                 .items():
            # output_files can be None if the connection is empty
            for output_file, input_files in output_files.items():
                if output_file != None and input_files != None:
                    max_level = max(
                        max_level, up_to_dateness_level(output_file))

        if max_level == 0:
            return self.get_step().get_pipeline().states.FINISHED
        elif max_level == 1:
            return self.get_step().get_pipeline().states.READY
        else:
            return self.get_step().get_pipeline().states.WAITING



    def add_private_info(self, key, value):
        '''
        Add private information to a run. Use this to store data which you will
        need when the run is executed. As opposed to public information,
        private information is not visible to subsequent steps.
        
        You can store paths to input files here, but not paths to output files as
        their expected location is not defined until we're in
        *AbstractStep.execute*
        (hint: they get written to a temporary directory inside *execute()*).
        '''
        if key in self._private_info and value != self._private_info[key]:
            logger.error(
                "You're trying to overwrite private info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._private_info[key]))
            sys.exit(1)
        self._private_info[key] = value

    def add_public_info(self, key, value):
        '''
        Add public information to a run. For example, a FASTQ reader may store 
        the index barcode here for subsequent steps to query via 
        ``AbstractStep.find_upstream_info()``.
        '''
        if key in self._public_info and value != self._public_info[key]:
            logger.error(
                "You're trying to overwrite public info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._public_info[key]))
            sys.exit(1)
        self._public_info[key] = value

    def update_public_info(self, key, value):
        '''
        Update public information already existing in a run. For example, all
        steps which handle FASTQ files want to know how to distinguish between
        files of read 1 and files of read 2. So each step that provides FASTQ
        should update this information if the file names are altered.
        The stored information can be acquired via:
        ``AbstractStep.find_upstream_info()``.
        '''
        if not key in self._public_info:
            logger.error("The key %s doesn't exist yet as public info."
                         "Please use add_public_info(%s, %s)"
                         % (key, key, value))
            sys.exit(1)
        else:
            self._public_info[key] = value

    def add_output_file(self, tag, out_path, in_paths):
        '''
        Add an output file to this run. Output file names must be unique across
        all runs defined by a step, so it may be a good idea to include the 
        run_id into the output filename.
        - *tag*: You must specify the connection annotation which must have been
                 previously declared via *AbstractStep.add_connection("out/...")*,
                 but this doesn't have to be done in the step constructor, it's
                 also possible in *declare_runs()* right before this method is
                 called.
        - *out_path*: The output file path, without a directory. The pipeline 
                      assigns directories for you (this parameter must not 
                      contain a slash).
        - *in_paths*: A list of input files this output file depends on. It is
                      **crucial** to get this right, so that the pipeline can
                      determine which steps are up-to-date at any given time.
                      You have to specify absolute paths here, including a 
                      directory, and you can obtain them via 
                      *AbstractStep.run_ids_and_input_files_for_connection*
                      and related functions.
        '''
        head, tail = os.path.split(out_path)

        # make sure there's no slash in out_path unless it's a source step
        if head != "" and not \
           isinstance(self._step, abst.AbstractSourceStep):
            logger.error("The declared output file path contains "
                         "directory separator: %s." % out_path)
            sys.exit(1)
        # make sure tag was declared with an outgoing connection
        if 'out/' + tag not in self._step._connections:
            logger.error("Invalid output_file tag '%s' in %s. "
                         "You might want to add self.add_connection('out/%s') "
                         "to the constructor of %s."
                         % (tag, str(self._step), tag, self._step.__module__))
            sys.exit(1)

        out_connection = self.get_out_connection(tag)

        if out_path in self.get_output_files_for_out_connection(out_connection):
            logger.error(
                "You're trying to re-add an output file which has already "
                "been declared: %s." % out_path)
            sys.exit(1)

        if not isinstance(in_paths, list):
            logger.error("Input paths (%s) is not a list." % in_paths)
            sys.exit(1)

        if None in in_paths:
            logger.error(
                "There is a NoneType element in input paths (%s) for output "
                "file (%s)" % (in_paths, out_path))
            sys.exit(1)

        if out_path == None:
            logger.error(
                "Trying to add NoneType element as output file for input paths "
                ": %s" % in_paths)
            sys.exit(1)

        self._input_files.append(in_paths)
        self._output_files[out_connection][out_path] = in_paths
        return_value = os.path.join(
                self.get_output_directory_du_jour_placeholder(), out_path)
        if head != "":
            return_value = os.path.abspath(out_path)
        return return_value

    @replace_output_dir_du_jour
    def add_temporary_file(self, prefix = '', suffix = '', designation = None):
        '''
        Returns the name of a temporary file (created by tempfile library).
        Name and output directory placeholder are concatenated. The concatenated
        string is returned and stored in a list. The placeholder is immediately
        properly adjusted by @replace_output_dir_du_jour.
        '''
        count = len(self._temp_paths)
        temp_placeholder = str()

        while True:
            hashtag = misc.str_to_sha1_b62('%s.%s.%s' % (prefix, count, suffix))
            temp_name = prefix + hashtag + suffix
            temp_placeholder = os.path.join(
                self.get_output_directory_du_jour_placeholder(), temp_name)

            if not temp_placeholder in self._temp_paths:
                break
            else:
                count += 1
        

        logger.info("Temporary file (#%s): %s" %
              (len(self._temp_paths) + 1, temp_name) )

        # _known_paths dict is logged
        known_paths = dict()
        known_paths[temp_placeholder] = {
            'label': os.path.basename(temp_placeholder),
            'designation': designation,
            'type': ''
        }
        self.add_known_paths(known_paths)
        # _temp_paths list contains all temporary files which are going to be
        # deleted
        self._temp_paths.append(temp_placeholder)
        return temp_placeholder

    def add_temporary_directory(self, prefix = '', suffix = '',
                                designation = None ):
        '''
        Convenience method for creation of temporary directories.
        Basically, just calls self.add_temporary_file().
        The magic happens in ProcessPool.__exit__()
        '''
        return self.add_temporary_file(prefix = prefix, suffix = suffix,
                                       designation = designation)

    def remove_temporary_paths(self):
        '''
        Everything stored in self._temp_paths is examined and deleted if
        possible. Also, self._known_paths 'type' info is updated here.
        NOTE: Included additional stat checks to detect FIFOs as well as other
        special files.
        '''
        for _ in self.get_temp_paths():
            # Check file type
            pathmode = os.stat(_).st_mode
            isdir = False if stat.S_ISDIR(pathmode) == 0 else True
            ischaracter = False if stat.S_ISCHR(pathmode) == 0 else True
            isblock = False if stat.S_ISBLK(pathmode) == 0 else True
            isfile = False if stat.S_ISREG(pathmode) == 0 else True
            isfifo = False if stat.S_ISFIFO(pathmode) == 0 else True
            islink = False if stat.S_ISLNK(pathmode) == 0 else True
            issock = False if stat.S_ISSOCK(pathmode) == 0 else True
            # Update 'type' value
            if _ in self.get_known_paths().keys():
                if isfile:
                    logger.debug("Set %s 'type' info to 'file'" % _)
                    self.get_known_paths()[_]['type'] = 'file'
                elif isdir:
                    logger.debug("Set %s 'type' info to 'directory'" % _)
                    self.get_known_paths()[_]['type'] = 'directory'
                elif isfifo:
                    logger.debug("Set %s 'type' info to 'fifo'" % _)
                    self.get_known_paths()[_]['type'] = 'fifo'
            if os.path.isdir(_) and isdir:
                try:
                    logger.info("Now deleting directory: %s" % _)
                    os.rmdir(_)
                except OSError as e:
                    logger.error("errno: %s" % e.errno)
                    logger.error("strerror: %s" % e.strerror)
                    logger.error("filename: %s" % e.filename)
                    pass
            else:
                try:
                    logger.info("Now deleting: %s" % _)
                    os.unlink(_)
                except OSError as e:
                    pass

    def add_empty_output_connection(self, tag):
        '''
        An empty output connection has 'None' as output file and 'None' as input
        file.
        '''
        # make sure tag was declared with an outgoing connection
        if 'out/' + tag not in self._step._connections:
            logger.error("Invalid output_file tag '%s' in %s. "
                         "You might want to add self.add_connection('out/%s') "
                         "to the constructor of %s."
                         % (tag, str(self._step), tag, self._step.__module__))
            sys.exit(1)
        try:
            out_connection = self.get_out_connection(tag)
        except KeyError:
            out_connection = self.add_out_connection(tag)

        if None in self._output_files[out_connection]:
            logger.error(
                "You're trying to re-declare %s as an empty output connection "
                % out_connection)
            sys.exit(1)

        self._output_files[out_connection][None] = None

    def add_out_connection(self, out_connection):
        if not out_connection.startswith("out/"):
            out_connection = 'out/' + out_connection
        self._output_files[out_connection] = dict()
        return out_connection

    def get_input_files_for_output_file(self, output_file):
        for connection in self.get_out_connections():
            if output_file in \
               self.get_output_files_for_out_connection(connection):
                return self._output_files[connection][output_file]

    def get_input_files_for_output_file_abspath(self, output_file):
        for connection in self.get_out_connections():
            if abspath_output_file in \
               self.get_output_files_abspath_for_out_connection(connection):
                return self.get_output_files_abspath()[connection]\
                    [abspath_output_file]

    def get_output_files_for_out_connection(self, out_connection):
        return list( self._output_files[out_connection].keys() )

    def get_output_files_abspath_for_out_connection(self, out_connection):
        return sorted(
            list( self.get_output_files_abspath()[out_connection].keys() )
        )

    def get_output_files(self):
        return self._output_files

    def get_output_files_abspath(self):
        '''
        Return a dictionary of all defined output files, grouped by connection 
        annotation::
        
           annotation_1:
               out_path_1: [in_path_1, in_path_2, ...]
               out_path_2: ...
           annotation_2: ...

        The ``out_path`` consists of the output directory du jour and the output
        file name.
        '''
        result = dict()
        for connection in self._output_files.keys():
            result[connection] = dict()
            for out_path, in_paths in self._output_files[connection].items():
                directory = self.get_output_directory_du_jour()
                full_path = out_path
                try:
                    head, tail = os.path.split(out_path)
                    if directory != None and out_path != None and head == "":
                        full_path = os.path.join(directory, out_path)
                except AttributeError:
                    pass
                result[connection][full_path] = in_paths

        return result

    def get_single_output_file_for_annotation(self, annotation):
        '''
        Retrieve exactly one output file of the given annotation, and crash
        if there isn't exactly one.
        '''
        temp = self.get_output_files_abspath()
        if len(temp[annotation]) != 1:
            logger.error("More than one output file declared for out/%s."
                         % annotation)
            sys.exit(1)
        return temp[annotation].keys()[0]

    def get_output_files_for_annotation_and_tags(self, annotation, tags):
        '''
        Retrieve a set of output files of the given annotation, assigned to
        the same number of specified tags. If you have two 'alignment' output
        files and they are called *out-a.txt* and *out-b.txt*, you can use this
        function like this:
        
        - *tags*: ['a', 'b']
        - result: {'a': 'out-a.txt', 'b': 'out-b.txt'}
        '''
        temp = self.get_output_files_abspath()
        return misc.assign_strings(temp[annotation].keys(), tags)

    def get_input_files_for_output_file(self, out_path):
        '''
        Return all input files a given output file depends on.
        '''
        temp = self.get_output_files()
        for tag in temp.keys():
            if out_path in temp[tag].keys():
                return sorted(temp[tag][out_path])
        logger.error("Sorry, your output '%s' file couldn't be found in"
                     "the dictionary: %s." % (out_path, temp))
        sys.exit(1)

    def get_public_info(self, key):
        '''
        Query public information which must have been previously stored via "
        "*add_public_info()*.
        '''
        return self._public_info[key]

    def has_public_info(self, key):
        '''
        Query whether a piece of public information has been defined.
        '''
        return (key in self._public_info)

    def get_private_info(self, key):
        '''
        Query private information which must have been previously stored via "
        "*add_private_info()*.
        '''
        return self._private_info[key]

    def has_private_info(self, key):
        '''
        Query whether a piece of public information has been defined.
        '''
        return (key in self._private_info)

    def as_dict(self):
        result = dict()
        result['output_directory'] = self.get_output_directory()
        result['output_files'] = self._output_files
        result['private_info'] = self._private_info
        result['public_info'] = self._public_info
        result['run_id'] = self._run_id
        return result

    def write_annotation_file(self, path):
        '''
        Write the YAML annotation after a successful or failed run. The
        annotation can later be used to render the process graph.
        '''
        
        # now write the annotation
        log = {}
        log['pid'] = os.getpid()
        log['step'] = {}
        log['step']['options'] = self.get_step().get_options()
        log['step']['name'] = self.get_step().get_step_name()
        log['step']['cores'] = self.get_step()._cores
        log['run'] = {}
        log['run']['run_info'] = self.as_dict()
        log['run']['run_id'] = self.get_run_id()
        log['run']['temp_directory'] = self.get_temp_output_directory()
        # if a submit script was used ...
        if os.path.exists(self.get_submit_script_file()):
            # ... read it and store it ...
            with open(self.get_submit_script_file(), 'r') as f:
                log['run']['submit_script'] = f.read()
            # ... finally delete it
            os.unlink(self.get_submit_script_file())
        log['run']['known_paths'] = self.get_known_paths()
        log['config'] = self.get_step().get_pipeline().config
        log['git_hash_tag'] = self.get_step().get_pipeline().git_hash_tag
        log['tool_versions'] = {}
        for tool in self.get_step()._tools.keys():
            log['tool_versions'][tool] = self.get_step().get_pipeline()\
                                                        .tool_versions[tool]
        log['pipeline_log'] = self.get_step()._pipeline_log
        log['start_time'] = self.get_step().start_time
        log['end_time'] = self.get_step().end_time
        if self.get_step().get_pipeline().git_dirty_diff:
            log['git_dirty_diff'] = self.get_step().get_pipeline().git_dirty_diff
        if self.get_step().get_pipeline().caught_signal is not None:
            log['signal'] = self.get_step().get_pipeline().caught_signal

        annotation_yaml = yaml.dump(log, default_flow_style = False)
        annotation_path = os.path.join(
            path, ".%s-annotation-%s.yaml" % 
            (self.get_run_id(), misc.str_to_sha1_b62(annotation_yaml)[:6]))

        # overwrite the annotation if it already exists
        with open(annotation_path, 'w') as f:
            f.write(annotation_yaml)
            
        return annotation_path, annotation_yaml
