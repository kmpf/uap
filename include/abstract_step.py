'''
Classes AbstractStep and AbstractSourceStep are defined here.

The class AbstractStep has to be inherited by all processing step classes.
The class AbstractSourceStep has to be inherited by all source step classes.

Processing steps generate output files from input files whereas source steps
only provide output files. Both step types may generates tasks, but only source
steps can introduce files from outside the destination path into the pipeline.
'''

# 1. standard library imports
import sys
import copy
import datetime
import hashlib
import inspect
from logging import getLogger
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
from shutil import copyfile
# 2. related third party imports
import fscache
import psutil
import yaml
# 3. local application/library specific imports
from uaperrors import UAPError
from connections_collector import ConnectionsCollector
import command as command_info
import misc
import process_pool
import pipeline_info
import run as run_module

abs_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(abs_path, 'steps'))
sys.path.insert(0, os.path.join(abs_path, 'sources'))
logger=getLogger('uap_logger')

class AbstractStep(object):

    fsc = fscache.FSCache()

    PING_TIMEOUT = 300
    PING_RENEW = 30
    VOLATILE_SUFFIX = '.volatile.placeholder.yaml'
    UNDERSCORE_OPTIONS = ['_depends', '_volatile', '_BREAK', '_connect',
                          '_cluster_submit_options', '_cluster_pre_job_command',
                          '_cluster_post_job_command', '_cluster_job_quota']

    states = misc.Enum(['DEFAULT', 'DECLARING', 'EXECUTING', 'POSTPROCESS'])

    def __init__(self, pipeline):

        self._pipeline = pipeline

        self.dependencies = list()
        '''
        All steps this step depends on.
        '''

        self._options = dict()
        '''
        Options as specified in the configuration.
        '''

        self._step_name = self.__module__
        '''
        By default, this is the name of the module. Can be overridden
        to allow for multiple steps of the same kind.
        '''

        self._runs = None
        '''
        Cached run information. ``declare_runs`` is only called once, the
        post-processed run objects are stored in here.
        '''

        self._pipeline_log = dict()

        self._cores = 1
        self._connections = set()
        self._optional_connections = set()
        self._connection_formats = dict()
        self._connection_descriptions = dict()
        self._connection_restrictions = dict()
        self._pre_command = dict()
        self._post_command = dict()
        self._module_load = dict()
        self._module_unload = dict()
        self._tools = dict()

        self._defined_options = dict()

        self.needs_parents = False

        self.children_step_names = set()

        self.finalized = False

        self._state = AbstractStep.states.DEFAULT

        self._submit_script = None

    def finalize(self):
        '''Finalizes the step.

        The intention is to make further changes to the step
        impossible, but apparently, it's checked nowhere at the moment.
        '''
        if self.finalized:
            return

        for parent_step in self.dependencies:
            parent_step.finalize()

        self.finalized = True

    def _reset(self):
        self._pipeline_log = dict()

    def get_pipeline(self):
        return self._pipeline

    def declare_run(self, run_id):
        '''
        Declare a run. Use it like this::

            with self.declare_run(run_id) as run:
                # add output files and information to the run here
        '''
        # Replace whitespaces by underscores
        run_id = re.sub(r'\s', '_', run_id)
        if run_id in self._runs:
            raise UAPError(
                "Cannot declare the same run ID twice: %s." % run_id)
        run = run_module.Run(self, run_id)
        self.add_run(run)
        return run

    def add_run(self, run):
        self._runs[run.get_run_id()] = run

    def get_run(self, run_id):
        '''
        Returns a single run object for run_id or None.
        '''
        if run_id in self._runs:
            return self._runs[run_id]
        else:
            return None

    def set_step_name(self, step_name):
        '''
        Change the step name.

        The step name is initially set to the module name. This method
        is used in case we need multiple steps of the same kind.
        '''
        self._step_name = step_name

    def set_options(self, options):
        '''
        Checks and stores step options.

        The options are either set to values given in YAML config or
        the default values set in self.add_option().
        '''
        self._options = dict()

        # set options
        for key, value in options.items():
            if key[0] == '_':
                if not key in AbstractStep.UNDERSCORE_OPTIONS:
                    raise UAPError(
                        "Invalid option in %s: %s" % (key, self))
                self._options[key] = value
            else:
                if not key in self._defined_options:
                    logger.info("Available options for %s:\n%s" %
                            (self.get_step_type(), yaml.dump(self._defined_options)))
                    raise UAPError(
                        "Unknown option in %s (%s): %s." %
                        (self.get_step_name(), self.get_step_type(), key))
                if value is not None and type(value) not in self._defined_options[key]['types']:
                    raise UAPError(
                        "Invalid type for option %s - it's %s and should be "
                        "one of %s." % (key, type(value),
                                        self._defined_options[key]['types']))
                if self._defined_options[key]['choices'] != None and \
                   value not in self._defined_options[key]['choices']:
                    raise UAPError(
                        "Invalid value '%s' specified for option %s - "
                        "possible values are %s." %
                        (value, key, self._defined_options[key]['choices']))
                self._options[key] = value

        # set default values for unset options and make sure all required
        # options have been set
        for key, info in self._defined_options.items():
            if key not in self._options:
                if info['optional'] is not True:
                    raise UAPError(
                        "Required option not set in %s: %s." % (self, key))
                self._options[key] = info['default']

        self._options.setdefault('_volatile', False)

        for i in ['_cluster_submit_options', '_cluster_pre_job_command',
                  '_cluster_post_job_command']:
            self._options.setdefault(i, '')
        self._options.setdefault('_cluster_job_quota', 0)

        self._options.setdefault('_connect', dict())
        self._options.setdefault('_depends', list())
        if not isinstance(self._options['_depends'], list):
            self._options['_depends'] = [self._options['_depends']]
        self._options['_depends'] = set(self._options['_depends'])
        # add implied dependencies
        for in_cons in self._options['_connect'].values():
            in_cons = in_cons if isinstance(in_cons, list) else [in_cons]
            for parent_cons in in_cons:
                parent = parent_cons.split("/")[0]
                self._options['_depends'].add(parent)

    def get_options(self):
        '''
        Returns a dictionary of all given options
        '''
        return self._options

    def get_option(self, key):
        """
        Query an option.
        """
        if key not in self._defined_options:
            raise UAPError(
                "Cannot query undefined option %s in step %s." %
                (key, self.__module__))
        return self._options[key]

    def is_option_set_in_config(self, key):
        """
        Determine whether an optional option (that is, a non-required option)
        has been set in the configuration.
        """
        if key not in self._defined_options:
            raise UAPError(
                "Cannot query undefined option %s in step %s." %
                (key, self.get_step_name()))
        is_set = key in self._options
        if is_set:
            if isinstance(self._options[key], list):
                is_set = any([v is not None for v in self._options[key]])
            else:
                is_set = self._options[key] is not None
        return is_set

    def is_volatile(self):
        return self._options['_volatile']

    def add_dependency(self, parent):
        '''
        Add a parent step to this steps dependencies.

        parent -- parent step this step depends on
        '''
        if not isinstance(parent, AbstractStep):
            raise UAPError(
                "Error: parent argument must be an AbstractStep.")
        if parent == self:
            raise UAPError("Cannot add a node as its own dependency.")
        self.dependencies.append(parent)
        parent.children_step_names.add(str(self))

    def get_dependencies(self):
        return self.dependencies

    def which_extensions_match_file_path(self, filepath, extensions):
        # Kann evtl. auch weg!
        extension_list = list()
        if type(filepath) is not str:
            raise UAPError("Filename must be string. Got %s of type %s"
                         % (filepath, type(filepath)))
        for ext in extensions:
            if type(ext) is not str:
                raise UAPError("Found non-string file extension: %s " % ext)
            else:
                extension_list.append(ext)

        ext_in_filename = list()
        file_parts = os.path.basename(filepath).split(".")
        for ext in extension_list:
            if ext in file_parts:
                ext_in_filename.append(ext)
        return ext_in_filename

    def get_input_runs(self):
        '''
        Return a dict which contains all runs per parent steps.
        '''
        input_runs = dict()
        for parent in self.get_dependencies():
            input_runs[parent.get_step_name()] = parent.get_runs()
        return input_runs

    def declare_runs(self):
        # fetch all incoming run IDs which produce reads...
        self.runs( self.get_run_ids_in_connections_input_files() )
        self.check_required_out_connections()

    def check_required_out_connections(self):
        '''
        This functions tests if all required out connections
        were set by all runs.
        '''
        required_out = self.get_out_connections(with_optional=False)
        bad_runs = 0
        for run_id, run in self._runs.items():
            used_conns = set()
            for connection, content in run._output_files.items():
                used = any([fl is not None for fl in content.keys()])
                if used:
                    used_conns.add(connection)
            missings = required_out - used_conns
            if missings:
                bad_runs += 1
                logger.warn('Run "%s" of step "%s" misses the required '
                        'connections %s. To remove this warning pass '
                        'optional=True to the add_connection method in the '
                        'step constructor __init__ of "%s".' %
                        (run_id, self.get_step_name(), list(missings),
                                self.get_step_type()))
            if bad_runs == 5:
                logger.warn('... Emitting connection test for further '
                    'runs of "%s".' % self.get_step_name())
                break
        if bad_runs:
            logger.warn('[Deprecation] Unmet required connections '
                    'may trigger an error in future version of the UAP.')

    def get_output_directory(self):
        '''
        Returns the step output directory.
        '''
        return os.path.join(
            self.get_pipeline().config['destination_path'],
            self.get_step_name()
        )

    def get_submit_script_file(self):
        if self._submit_script == None:
            self._submit_script = os.path.join(
                self.get_output_directory(),
                ".submit-%s.sh" % self.get_step_name()
            )
        return self._submit_script

    def runs(self, run_ids_connections_files):
        '''
        Abstract method this must be implemented by actual step.

        Raise NotImplementedError if subclass does not override this
        method.
        '''
        raise NotImplementedError()

    def execute(self, run_id, run):
        # get run_info objects
        with self.get_run(run_id) as run:
            logger.info("Run ID: %s" % run_id)
            # for each exec_group in that run ...
            for exec_group in run.get_exec_groups():
                # ... create a process pool
                with process_pool.ProcessPool(run) as pool:
                    # Clean up (use last ProcessPool for that)
                    if exec_group == run.get_exec_groups()[-1]:
                        logger.info("Telling pipeline to clean up!")
                        pool.clean_up_temp_paths()

                    for poc in exec_group.get_pipes_and_commands():
                        # for each pipe or command (poc)
                        # check if it is a pipeline ...
                        if isinstance(poc, pipeline_info.PipelineInfo):
                            # ... create a pipeline ...
                            with pool.Pipeline(pool) as pipeline:
                                for command in poc.get_commands():
                                    pipeline.append(
                                        command.get_command(),
                                        stdout_path = command.get_stdout_path(),
                                        stderr_path = command.get_stderr_path())
                        elif isinstance(poc, command_info.CommandInfo):
                            pool.launch(
                                poc.get_command(),
                                stdout_path = poc.get_stdout_path(),
                                stderr_path = poc.get_stderr_path())

    def get_runs(self):
        '''
        Getter method for runs of this step.

        If there are no runs as this method is called, they are created here.
        '''
        # create runs if they don't exist yet
        if not self._runs:
            # if _BREAK: true is specified in the configuration,
            # return no runs and thus cut off further processing
            if '_BREAK' in self._options and self._options['_BREAK']:
                return dict()

            self._runs = dict()

            self._state = AbstractStep.states.DECLARING
            self.declare_runs()
            self._state = AbstractStep.states.DEFAULT

            # define file dependencies
            for run_id in self._runs.keys():
                pipeline = self.get_pipeline()
                run = self.get_run(run_id)
                for connection in run.get_output_files_abspath().keys():
                    for output_path, input_paths in \
                        run.get_output_files_abspath()[connection].items():
                        # proceed if we have normal output_path/input_paths
                        if output_path != None and input_paths != None:
                            # store file dependencies
                            pipeline.add_file_dependencies(
                                output_path, input_paths)
                            # create task ID
                            task_id = '%s/%s' % (str(self), run_id)
                            pipeline.add_task_for_output_file(
                                output_path, task_id)
                            # No input paths? Add empty string NOT None
                            # as file name
                            if len(input_paths) == 0:
                                pipeline.add_task_for_input_file(
                                    "", task_id)
                            for input_path in input_paths:
                                pipeline.add_task_for_input_file(
                                    input_path, task_id)

        # now that _runs exists, it remains constant, just return it
        return self._runs

    def get_run_ids(self):
        '''
        Returns sorted list of runs generated by step.
        '''
        return sorted(self.get_runs().keys())

    def get_step_name(self):
        '''
        Returns this steps name.

        Returns the step name which is initially equal to the step type
        (== module name)  but can be changed via set_step_name() or via
        the YAML configuration.
        '''
        return self._step_name

    def get_step_type(self):
        '''
        Returns the original step name (== module name).
        '''
        return self.__module__

    def get_run_state_basic(self, run_id):
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
            It then checks all upstream paths, which may in turn be volatile
            placeholder files.
            '''

            # reconstruct original path from volatile placeholder path
            path = volatile_path[:-len(AbstractStep.VOLATILE_SUFFIX)]

            if AbstractStep.fsc.exists(path):
                # the original file still exists, ignore volatile placeholder
                return False

            if not path in self.get_pipeline().task_id_for_output_file:
                # there is no task which creates the output file
                return False

            task_id = self.get_pipeline().task_id_for_output_file[path]

            task = self.get_pipeline().task_for_task_id[task_id]
            if not task.step.is_volatile():
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
            if path in self.get_pipeline().file_dependencies_reverse:
                uncovered_files = self.get_pipeline().\
                                  file_dependencies_reverse[path]

            for downstream_path, downstream_info in info['downstream'].items():
                if downstream_path in self.get_pipeline().task_id_for_output_file:
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
                    if pv_downstream_path.endswith(AbstractStep.VOLATILE_SUFFIX):
                        if not volatile_path_good(pv_downstream_path, recurse):
                            return False

            if len(uncovered_files) > 0:
                # there are still files defined which are not covered by the
                # placeholder
                return False

            return True

        def change_to_volatile_if_need_be(path, recurse = True):
            """
            Changes the file path to volatile path if necessary."""
            if path != None:
                if not AbstractStep.fsc.exists(path):
                    # the real output file does not exist
                    volatile_path = path + AbstractStep.VOLATILE_SUFFIX
                    if volatile_path_good(volatile_path, recurse):
                        return volatile_path
                return path

        def is_path_up_to_date(outpath, inpaths):
            """
            First, replace paths with volatile paths if the step is marked
            as volatile and the real path is missing.
            But, only consider volatile placeholders if all child tasks
            are finished. That means if a child of a volatile
            step needs to be run because it has been added or an existing step
            has been modified, the volatile placeholders are ignored, thus
            turning the task from 'finished' to 'ready' or 'waiting'
            Hint: The pv_ prefix is for 'possibly volatile'
            """
            pv_outpath = outpath
            pv_inpaths = list()

            if outpath in self.get_pipeline().task_id_for_output_file:
                pv_outpath = change_to_volatile_if_need_be(outpath)

            for inpath in inpaths:
                pv_inpaths.append(change_to_volatile_if_need_be(inpath))

            if not AbstractStep.fsc.exists(pv_outpath):
                if not inpaths and AbstractStep.fsc.islink(pv_outpath):
                    # these are probably files of a broken source step
                    raise UAPError('The input "%s" seems to be a broken '
                                   'symlink.' % pv_outpath)
                return False
            for pv_inpath in pv_inpaths:
                if not AbstractStep.fsc.exists(pv_inpath):
                    return False
                # Check that inpath was last modified before outpath
                if AbstractStep.fsc.getmtime(pv_inpath) > \
                   AbstractStep.fsc.getmtime(pv_outpath):
                    logger.info('"%s" was changed befor its dependency "%s".' %
                            (pv_outpath, pv_inpaths))
                    return False
            return True

        def up_to_dateness_level(path, level = 0):
            result = level

#            sys.stderr.write("path: %s (level: %d)\n" % (path, level))

            if path != None:
                dep_paths = self.get_pipeline().file_dependencies[path]
#                sys.stderr.write("dep_paths: %s\n" % dep_paths)
                if not is_path_up_to_date(path, dep_paths):
                    result = level + 1
                for dep_path in dep_paths:
                    recursive_result = up_to_dateness_level(dep_path, level + 1)
                    if recursive_result > level + 1:
                        result = max(result, recursive_result)
                return result

        """
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
        """

        run_info = self.get_runs()
        max_level = 0
        for tag, output_files in run_info[run_id].get_output_files_abspath().items():

            # output_files can be None if the connection is empty
            for output_file, input_files in output_files.items():
                if output_file != None and input_files != None:
#                    sys.stderr.write("outputfile: %s\n" % output_file)
                    max_level = max(
                        max_level, up_to_dateness_level(output_file))

        if max_level == 0:
            return self.get_pipeline().states.FINISHED
        elif max_level == 1:
            return self.get_pipeline().states.READY
        else:
            return self.get_pipeline().states.WAITING

    def get_run_state(self, run_id):
        '''
        Returns run state of a run.

        Determine the run state (that is, not *basic* but *extended* run state)
        of a run, building on the value returned by *get_run_state_basic()*.

        If a run is **ready**, this will:
          - return **executing** if an up-to-date *executing ping file* is found
          - otherwise return **queued** if a *queued ping file* is found

        If a run is **waiting**, this will:
          - return **queued** if a *queued ping file* is found

        Otherwise, it will just return the value obtained from
        *get_run_state_basic()*.

        Attention: The status indicators **executing** and **queued** may be
        temporarily wrong due to the possiblity of having out-of-date ping files
        lying around.
        '''
        run = self.get_run(run_id)
        run_state = self.get_run_state_basic(run_id)
        if run_state == self.get_pipeline().states.READY:
            if AbstractStep.fsc.exists( run.get_executing_ping_file() ):
                # here, we just check whether the executing ping file exists,
                # it doesn't matter whether it's been stale for a year
                # (the user will get notified that there are stale ping files
                # and can fix it with ./fix-problems.py, it's probably better
                # to fix this explicitly
                return self.get_pipeline().states.EXECUTING
            if AbstractStep.fsc.exists( run.get_queued_ping_file() ):
                return self.get_pipeline().states.QUEUED
        elif run_state == self.get_pipeline().states.WAITING:
            if AbstractStep.fsc.exists( run.get_queued_ping_file() ):
                return self.get_pipeline().states.QUEUED
        return run_state

    def _move_ping_files(self, executing_ping_path, queued_ping_path):
        # don't remove the ping file, rename it so we can inspect it later
        suffix = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        if os.path.exists(executing_ping_path):
            try:
                out_w_suffix = executing_ping_path + '.' + suffix
                os.rename(executing_ping_path, out_w_suffix)
                logger.debug('The run ping file "%s" was moved to "%s" '
                             'by host %s.' %
                             (executing_ping_path, out_w_suffix,
                                     socket.gethostname()))
            except OSError as e:
                logger.debug('The run ping file "%s" could not be moved: %s' %
                             (executing_ping_path, str(e)))
                pass
        else:
            logger.debug('This run ping file was not found: %s' %
                         executing_ping_path)
        # remove the queued ping file
        if os.path.exists(queued_ping_path):
            try:
                out_w_suffix = queued_ping_path + '.' + suffix
                os.rename(queued_ping_path,
                          queued_ping_path + '.last')
                copyfile(queued_ping_path + '.last', out_w_suffix)
                logger.debug('The queued ping file "%s" recived the ".last" '
                             'suffix and was copied to "%s" by host %s.' %
                             (queued_ping_path, out_w_suffix,
                                     socket.gethostname()))
            except OSError as e:
                logger.debug('The queued ping file "%s" could not be moved '
                             'or copied: %s' % (queued_ping_path, str(e)))
                pass
        else:
            logger.debug('This queued ping file was not found: %s' %
                    queued_ping_path)


    def run(self, run_id):
        '''
        Create a temporary output directory and execute a run. After the run
        has finished, it is checked that all output files are in place and
        the output files are moved to the final output location. Finally,
        YAML annotations are written.
        '''

        # this is the run we'll execute now
        run = self.get_run(run_id)

        # create the output directory if it doesn't exist yet
        if not os.path.isdir(run.get_output_directory()):
            os.makedirs(run.get_output_directory())

        # now write the run ping file
        executing_ping_path = run.get_executing_ping_file()

        if os.path.exists(executing_ping_path):
            raise UAPError("%s/%s seems to be already running, exiting..."
                         % (self, run_id))
        queued_ping_path = run.get_queued_ping_file()

        # create a temporary directory for the output files
        temp_directory = run.get_temp_output_directory()
        os.makedirs(temp_directory)

        # prepare known_paths
        known_paths = dict()
        for tag, tag_info in run.get_output_files_abspath().items():
            for output_path, input_paths in tag_info.items():
                # add the real output path
                if output_path != None and input_paths != None:
                    known_paths[output_path] = {
                        'type': 'output',
                        'designation': 'output',
                        'label': os.path.basename(output_path),
                        'type': 'step_file'}
                    # ...and also add the temporary output path
                    known_paths[
                        os.path.join(temp_directory, os.path.basename(
                            output_path))] = {
                        'type': 'output',
                        'designation': 'output',
                        'label': "%s\\n(%s)" %
                                (os.path.basename(output_path), tag),
                        'type': 'step_file',
                        'real_path': output_path}
                    for input_path in input_paths:
                        if input_path != None:
                            known_paths[input_path] = {
                                'type': 'input',
                                'designation': 'input',
                                'label': os.path.basename(input_path),
                                'type': 'step_file'}

        # now write the run ping file
        executing_ping_info = dict()
        executing_ping_info['start_time'] = datetime.datetime.now()
        executing_ping_info['host'] = socket.gethostname()
        executing_ping_info['pid'] = os.getpid()
        executing_ping_info['cwd'] = os.getcwd()
        executing_ping_info['temp_directory'] = run.get_temp_output_directory()

        with open(executing_ping_path, 'w') as f:
            f.write(yaml.dump(executing_ping_info, default_flow_style = False))

        executing_ping_pid = os.fork()
        if executing_ping_pid == 0:
            # this is the chid process
            try:
                signal.signal(signal.SIGTERM, signal.SIG_DFL)
                signal.signal(signal.SIGINT, signal.SIG_IGN)
                while True:
                    time.sleep(AbstractStep.PING_RENEW)
                    # if the executing ping file is gone and the touching
                    # operation fails, then SO BE IT!
                    os.utime(executing_ping_path, None)
            finally:
                os._exit(0)

        original_term_handler = signal.getsignal(signal.SIGTERM)
        original_int_handler = signal.getsignal(signal.SIGINT)
        def ping_on_term(signum, frame):
            logger.debug('Recived SIGTERM and moving ping files...')
            self._move_ping_files(executing_ping_path, queued_ping_path)
            original_term_handler(signum, frame)
        def ping_on_int(signum, frame):
            logger.debug('Recived SIGINT and moving ping files...')
            self._move_ping_files(executing_ping_path, queued_ping_path)
            original_int_handler(signum, frame)

        signal.signal(signal.SIGTERM, ping_on_term)
        signal.signal(signal.SIGINT, ping_on_int)

        self.start_time = datetime.datetime.now()
        self.get_pipeline().notify(
            "[START] [%s] starting %s/%s on %s" %
            (self.get_pipeline().config['id'], str(self), run_id,
             socket.gethostname()))
        caught_exception = None
        self._state = AbstractStep.states.EXECUTING
        base_working_dir = os.getcwd()
        os.chdir(run.get_temp_output_directory())
        try:
            self.execute(run_id, run)
        except Exception as e:
            if type(e) != UAPError:
                logger.error("%s: %s" % (type(e).__name__, e))
            # Oh my. We have a situation. This is awkward. Tell the process
            # pool to wrap up. This way, we can try to get process stats before
            # shutting everything down.
            process_pool.ProcessPool.kill()
            # Store the exception, re-raise it later
            caught_exception = sys.exc_info()
        finally:
            signal.signal(signal.SIGTERM, original_term_handler)
            signal.signal(signal.SIGINT, original_int_handler)
            self._move_ping_files(executing_ping_path, queued_ping_path)
            self._state = AbstractStep.states.POSTPROCESS # changes relative paths
            os.chdir(base_working_dir)
            try:
                os.kill(executing_ping_pid, signal.SIGTERM)
                os.waitpid(executing_ping_pid, 0)
            except OSError:
                # if the ping process was already killed, it's gone anyway
                pass

        # TODO: Clean this up. Re-think exceptions and task state transisitions.

        self.end_time = datetime.datetime.now()

        if (not self.get_pipeline().caught_signal) and (caught_exception is None):
            # if we're here, we can assume the step has finished successfully
            # now rename the output files (move from temp directory to
            # destination directory)

            # import pdb
            # pdb.set_trace()

            for tag in run.get_output_files().keys():
                for out_file in run.get_output_files()[tag].keys():
                    # don't try to rename files if they were not meant to exist
                    # in our temporary directory
                    # 1. out_file should not be None (empty output connection)
                    # 2. out_file should not contain a '/' (file belongs to a
                    #    source step)
                    if out_file != None and not '/' in out_file:
                        source_path = os.path.join(
                            run.get_temp_output_directory(),
                            os.path.basename(out_file)
                        )
                        destination_path = os.path.join(
                            run.get_output_directory(),
                            os.path.basename(out_file))
                        # first, delete a possibly existing volatile placeholder
                        # file
                        destination_path_volatile = destination_path + \
                                                    AbstractStep.VOLATILE_SUFFIX
                        if os.path.exists(destination_path_volatile):
                            logger.info("Now deleting: %s" % destination_path_volatile)
                            os.unlink(destination_path_volatile)
                        # TODO: if the destination path already exists, this
                        # will overwrite the file.
                        if os.path.exists(source_path):
                            # Calculate SHA256 hash for output files
                            sha256sum = hashlib.sha256()
                            try:
                                with open(source_path, 'rb') as f:
                                    # the below exception is raised for large files
                                    # this workaround reads the file in chunks and
                                    # updates the sha256sum
                                    while True:
                                        # read file in 2MB chunks
                                        buf = f.read(2*1024*1024)
                                        if not buf:
                                            break
                                        sha256sum.update(buf)
                            except:
                                raise UAPError("Error while calculating SHA256sum "
                                             "of %s" % source_path)

                            # hexadecimal version of sha256sum
                            sha256sum = sha256sum.hexdigest()

                            os.rename(source_path, destination_path)
                            for path in [source_path, destination_path]:
                                if path in known_paths.keys():
                                    if known_paths[path]['designation'] == \
                                       'output':
                                       known_paths[path]['sha256'] = sha256sum
                                    if known_paths[path]['type'] != \
                                       'step_file':
                                        logger.debug("Set %s 'type' info to "
                                                     "'step_file'" % path)
                                        known_paths[path]['type'] = 'step_file'
                        else:
                            caught_exception = (
                                None,
                                StandardError(
                                    "The step failed to produce an announced "\
                                    "output file: %s. "\
                                    "\nSource file doesn't exists: %s" %
                                    (out_file, source_path)),
                                None)

        for path, path_info in known_paths.items():
            # Get the file size
            if os.path.exists(path):
                known_paths[path]['size'] = os.path.getsize(path)

        run.add_known_paths(known_paths)
        annotation_path, annotation_str = run.write_annotation_file(
            run.get_output_directory() \
            if ((self.get_pipeline().caught_signal is None) and \
                (caught_exception is None)) \
            else run.get_temp_output_directory())

        self._state = AbstractStep.states.DEFAULT

        if self.get_pipeline().caught_signal is not None or \
           caught_exception is not None:
            message = "[BAD] %s/%s failed on %s after %s\n" % \
                      (str(self), run_id, socket.gethostname(),
                       misc.duration_to_str(self.end_time - self.start_time))
            message += "Here are the details: " + annotation_path + '\n'
            attachment = None
            if os.path.exists(annotation_path + '.png'):
                attachment = dict()
                attachment['name'] = 'details.png'
                attachment['data'] = open(annotation_path + '.png').read()
            # todo as: remove the following line?
            self.get_pipeline().notify(message, attachment)
            if caught_exception is not None:
                raise caught_exception[1], None, caught_exception[2]
        else:
            # create a symbolic link to the annotation for every output file
            for tag in run._output_files.keys():
                for out_path in run._output_files[tag].keys():
                    if out_path != None:
                        destination_path = os.path.join(
                            run.get_output_directory(),
                            '.' + os.path.basename(out_path) +
                            '.annotation.yaml')
                        # overwrite the symbolic link if it already exists
                        if os.path.exists(destination_path):
                            logger.info("Now deleting: %s" % destination_path)
                            os.unlink(destination_path)
                        oldwd = os.getcwd()
                        os.chdir(os.path.dirname(destination_path))
                        os.symlink(os.path.basename(annotation_path),
                                   os.path.basename(destination_path))
                        os.chdir(oldwd)

            # finally, remove the temporary directory if it's empty
            try:
                os.rmdir(temp_directory)
            except OSError as e:
                logger.info('Coult not remove temp dir "%s": %s' %
                        (temp_directory, e))

            # step has completed successfully, now determine how many jobs are
            # still left but first invalidate the FS cache because things have
            # changed by now...
            AbstractStep.fsc = fscache.FSCache()

            remaining_task_info = self.get_run_info_str()

            message = "[OK] %s/%s successfully finished on %s after %s\n" % \
                      (str(self), run_id, socket.gethostname(),
                       misc.duration_to_str(self.end_time - self.start_time))
            message += str(self) + ': ' + remaining_task_info + "\n"
            attachment = None
            if os.path.exists(annotation_path + '.png'):
                attachment = dict()
                attachment['name'] = 'details.png'
                attachment['data'] = open(annotation_path + '.png').read()
            self.get_pipeline().notify(message, attachment)

            self._reset()

    def get_parents(self, run_id):
        '''
        Returns all tasks that produce input for this task.
        '''
        task_id = '%s/%s' % (self, run_id)
        input_files = set()
        if task_id in self.get_pipeline().input_files_for_task_id:
            input_files = self.get_pipeline().input_files_for_task_id[task_id]
        parents = set()
        # Only source steps do have empty strings in the input files list
        # so we can safely exclude them here
        for inpath in [x for x in input_files if x != '']:
            task_id = self.get_pipeline().task_id_for_output_file[inpath]
            if task_id in self.get_pipeline().task_for_task_id:
                task = self.get_pipeline().task_for_task_id[task_id]
                if task.step._options['_volatile'] == True:
                    parents.add(task)
        return parents

    def reports(self, run_id, out_connection_output_files):
        '''
        Abstract method this must be implemented by actual step.

        Raise NotImplementedError if subclass does not override this
        method.
        '''
        raise NotImplementedError()


    def generate_report(self, run_id):
        '''
        Gathers the output files for each outgoing connection and calls
        self.reports() to do the job of creating a report.
        '''

        run = self.get_run(run_id)
        out_connection_output_files = dict()
        for out_connection in run.get_out_connections():
            out_connection_output_files[out_connection] = run.\
                get_output_files_abspath_for_out_connection(out_connection)

        try:
            self.reports(run_id, out_connection_output_files)
        except NotImplementedError as e:
            logger.info('Step %s is not capable to generate reports' %
                        (self._step_name))
        except Exception as e:
            logger.error('Unexpected error while trying to generate report for '
                         'task %s/%s: %s' % (self._step_name, run_id,
                                             e))

    def generate_one_report(self):
        '''
        Gathers the output files for each outgoing connection and calls
        self.reports() to do the job of creating a report.
        '''

        run_ids_connections_output_files = self\
            .get_run_ids_out_connections_output_files()
        for run_id in run_ids_connections_output_files:
            for con, files in run_ids_connections_output_files[run_id].items():
                files = [f for f in files if os.path.isfile(f)]
                run_ids_connections_output_files[run_id][con] = files
        try:
            self.reports( run_ids_connections_output_files )
        except NotImplementedError as e:
            logger.info('Step %s is not capable to generate reports' %
                        (self._step_name))
        except Exception as e:
            raise UAPError('Unexpected error while trying to generate report for '
                         'step %s: %s' % (self.get_step_name(), e))

    def get_pre_commands(self):
        """
        Return dictionary with commands to execute before starting any other
        command of this step
        """
        return self._pre_command

    def get_module_loads(self):
        """
        Return dictionary with module load commands to execute before starting
        any other command of this step
        """
        return self._module_load

    def get_tool(self, key):
        """
        Return full path to a configured tool.
        """
        if key not in self._tools:
            raise UAPError("Tool %s unknown. Maybe you forgot to use "
                         "self.require_tool('%s')" % (key, key))
        return self._tools[key]

    def get_module_unloads(self):
        """
        Return dictionary with module unload commands to execute before
        starting any other command of this step
        """
        return self._module_unload


    def get_post_commands(self):
        """
        Return dictionary with commands to execute after finishing any other
        command of this step
        """
        return self._post_command


    def get_run_info_str(self):
        count = {}
        for _ in self.get_run_ids():
#            sys.stderr.write("cur_run_id: %s\n" % _)
            state = self.get_run_state(_)
#            sys.stderr.write("state: %s\n" % state)
            if not state in count:
                count[state] = 0
            count[state] += 1
        return ', '.join(["%d %s" % (count[_], _.lower()) \
                          for _ in self.get_pipeline().states.order if _ in count])

    def append_pipeline_log(self, log):
        if len(self._pipeline_log) == 0:
            self._pipeline_log = log
        else:
            for k in log.keys():
                if k == 'process_watcher':
                    for k2 in log[k].keys():
                        if k2 == 'max':
                            for _ in log[k][k2].keys():
                                if _ == 'sum':
                                    for k3 in self._pipeline_log[k][k2][_].keys():
                                        self._pipeline_log[k][k2][_][k3] = \
                                            max(self._pipeline_log[k][k2][_][k3],
                                                log[k][k2][_][k3])
                                else:
                                    self._pipeline_log[k][k2][_] = log[k][k2][_]
                        else:
                            self._pipeline_log[k][k2].update(log[k][k2])

                else:
                    if log[k].__class__ == list:
                        self._pipeline_log[k].extend(log[k])
                    else:
                        self._pipeline_log[k].update(log[k])


    def __str__(self):
        return self._step_name

    @classmethod
    def get_step_class_for_key(cls, key):
        """
        Returns a step (or source step) class for a given key which corresponds
        to the name of the module the class is defined in. Pass 'cutadapt' and
        you will get the cutadapt.Cutadapt class which you may then instantiate.
        """

        # Attention, import statement in class method coming right up!
        # Ok, this is strange, I know. But we need the io_step.IOStep class now
        # because we want to test whether module members are a subclass of this
        # and if we import it right at the beginning of this file, we would create
        # a circular reference, because AbstractStep is imported at the beginning
        # of io_step. There's probably a better solution, but I think it doesn't
        # hurt, either. Here goes the awkward line:


        check_classes = [AbstractSourceStep, AbstractStep]
        for index, c in enumerate(check_classes):

            classes = [_ for _ in inspect.getmembers(__import__(key),
                                                     inspect.isclass) \
                       if c in _[1].__bases__]

            for k in range(index):
                classes = [_ for _ in classes if _[1] != check_classes[k]]
            if len(classes) > 0:
                if len(classes) != 1:
                    raise UAPError("need exactly one subclass of %s in %s"
                                 % (c, key))
                return classes[0][1]

        raise UAPError("No suitable class found for module %s." % key)
    def set_cores(self, cores):
        """
        Specify the number of CPU cores this step will use.
        """
        self._cores = cores

    def get_cores(self):
        """
        Returns the number of cores used in this step.
        """
        return self._cores

    def add_input_connection(self, connection, constraints = None):
        '''
        Add an input connection to this step
        '''
        self.add_connection('in/%s' % connection, constraints)

    def add_output_connection(self, connection, constraints = None):
        '''
        Add an output connection to this step
        '''
        self.add_connection('out/%s' % connection, constraints)

    def add_connection(self, connection, constraints = None,
            optional = False, format = None, description = None):
        """
        Add a connection, which must start with 'in/' or 'out/'.
        :type format: (str) Data format passed in the connection.
        :type description: (str) Explain the connection.
        """
        if not (connection[0:3] == 'in/' or connection[0:4] == 'out/'):
            raise UAPError("A connection must start with 'in/' or 'out/'.")
        if connection[0:3] == 'in/':
            self.needs_parents = True
        if optional is True:
            self._optional_connections.add(connection)
        else:
            self._connections.add(connection)
        if format is not None:
            self._connection_formats[connection] = format
        if description is not None:
            self._connection_descriptions[connection] = description
        if constraints is not None:
            self._connection_restrictions[connection] = constraints

    def get_connections(self, with_optional=True):
        """
        Return all connections for this step
        """
        connections = self._connections.copy()
        if with_optional is True:
            connections = connections.union(self._optional_connections)
        return connections

    def get_in_connections(self, with_optional=True, strip_prefix=False):
        """
        Return all in-connections for this step
        """
        connections = self._connections.copy()
        if with_optional is True:
            connections = connections.union(self._optional_connections)
        in_connections = set()
        for connection in connections:
            if connection[0:3] == "in/":
                if strip_prefix is True:
                    con = connection[3:]
                else:
                    con = connection
                in_connections.add(con)
        return in_connections

    def get_out_connections(self, with_optional=True, strip_prefix=False):
        """
        Return all out-connections for this step
        """
        connections = self._connections.copy()
        if with_optional is True:
            connections = connections.union(self._optional_connections)
        out_connections = set()
        for connection in connections:
            if connection[0:4] == "out/":
                if strip_prefix is True:
                    con = connection[4:]
                else:
                    con = connection
                out_connections.add(con)
        return out_connections

    def require_tool(self, tool):
        """
        Declare that this step requires an external tool. Query it later with
        *get_tool()*.
        """
        if self.get_pipeline() is not None:
            if not tool in self.get_pipeline().config['tools']:
                raise UAPError("%s requires the tool %s but it's not declared in "
                             "the configuration."
                             % (self, tool))
            self._tools[tool] = self.get_pipeline().config['tools'][tool]['path']
            if 'pre_command' in self.get_pipeline().config['tools'][tool]:
                self._pre_command[tool] = self.get_pipeline().config['tools'][tool]\
                                          ['pre_command']
            if 'module_load' in self.get_pipeline().config['tools'][tool]:
                self._module_load[tool] = self.get_pipeline().config['tools'][tool]\
                                          ['module_load']
            if 'module_unload' in self.get_pipeline().config['tools'][tool]:
                self._module_unload[tool] = self.get_pipeline().config['tools'][tool]\
                                            ['module_unload']
            if 'post_command' in self.get_pipeline().config['tools'][tool]:
                self._post_command[tool] = self.get_pipeline().config['tools'][tool]\
                                           ['post_command']
        else:
            self._tools[tool] = True

    def add_option(self, key, *option_types, **kwargs):
        """
        Add an option. Multiple types may be specified.
        """
        if not 'optional' in kwargs:
            kwargs['optional'] = False
        for _ in ['default', 'description', 'choices']:
            if not _ in kwargs:
                kwargs[_] = None

        if key[0] == '_':
            raise UAPError(
                "Option key must not start with an underscore: %s." % key)
        if key in self._defined_options:
            raise UAPError("Option %s is already defined." % key)
        if len(option_types) == 0:
            raise UAPError("No option type specified for option %s." % key)
        if len(option_types) > 1 and kwargs['choices'] != None:
            raise UAPError("You cannot define choices if multiple options types "
                         "are defined (%s)." % key)
        for option_type in option_types:
            if not  option_type in [int, float, str, bool, list, dict]:
                raise UAPError("Invalid type for option %s: %s."
                             % (key, option_type))
        if kwargs['optional'] and (kwargs['default'] != None):
            if type(kwargs['default']) not in option_types:
                raise UAPError(
                    "In step: (%s) option: (%s) Type of default value (%s) does not match any of the "
                    "declared possible types (%s)."
                    % (self, key, type(kwargs['default']), option_types))

        info = dict()
        info['types'] = option_types
        for _ in ['optional', 'default', 'description', 'choices']:
            info[_] = kwargs[_]

        if info['description'] is not None:
            if not isinstance(info['description'], str):
                raise UAPError('The description of option %s in step %s is not a string.' %
                        (key, self))
            # collapse whites spaces
            info['description'] = re.sub('\s+', ' ', info['description'])

        self._defined_options[key] = info

    def find_upstream_info_for_input_paths_as_set(self, input_paths,
                                                  key, expected = 1):
        task_ids = set()
        for path in input_paths:
            task_ids.add(self.get_pipeline().task_id_for_output_file[path])
        results = set()
        for task_id in task_ids:
            task = self.get_pipeline().task_for_task_id[task_id]
            step = task.step
            run_id = task.run_id
            run = step._runs[run_id]
            if run.has_public_info(key):
            	results.add(run.get_public_info(key))
            results |= self.find_upstream_info_for_input_paths_as_set(
                task.input_files(), key, None)

        if expected is not None:
            if len(results) != expected:
                raise UAPError(
                    "Unable to determine upstream %s info from %s." %
                    (key, self))
        return results

    def find_upstream_info_for_input_paths(self, input_paths, key):
        """
        Find a piece of public information in all upstream steps. If the
        information is not found or defined in more than one upstream step,
        this will crash.
        """
        # And boy, will it crash. SUH-MAAAASH! http://youtu.be/PbYD7sj6vxc?t=1m38s

        result = self.find_upstream_info_for_input_paths_as_set(
            input_paths, key, expected = 1)
        return list(result)[0]

    def get_run_ids_out_connections_output_files(self):
        '''
        Return a dictionary with all run IDs of the current step, their
        out connections, and the files that belong to them::

           run_id_1:
               in_connection_1: [input_path_1, input_path_2, ...]
               in_connection_2: ...
           run_id_2: ...

        Format of ``in_connection``: ``in/<connection>``. Input paths are
        absolute.
        '''
        run_ids_connections_files = dict()

        for run in self.get_runs():
            run_id = run.get_run_id()
            run_ids_connections_files[run_id] = dict()
            for out_connection in run.get_out_connections():
                run_ids_connections_files[run_id][out_connection] = run\
                    .get_output_files_for_out_connection(out_connection)

        return run_ids_connections_files


    def get_run_ids_in_connections_input_files(self):
        '''
        Return a dictionary with all run IDs from parent steps, the
        in connections they provide data for, and the names of the files::

           run_id_1:
               in_connection_1: [input_path_1, input_path_2, ...]
               in_connection_2: ...
           run_id_2: ...

        Format of ``in_connection``: ``in/<connection>``. Input paths are
        absolute.
        '''

        cc = ConnectionsCollector(self.get_step_name())
        self._options.setdefault('_connect', dict())

        # Check if set in-connections are defined in the step class
        # and collect out connections for later check.
        set_out_connections = set()
        used_out_connections = set()
        for in_conn, out_conn in self._options['_connect'].items():
            if in_conn not in self.get_in_connections():
                raise UAPError('_connect: unknown input connection "%s" '
                             'found.' % in_conn)
            out_conn = out_conn if isinstance(out_conn, list) else [out_conn]
            set_out_connections = set_out_connections.union(set(out_conn))

        if 'empty' in set_out_connections:
            logger.warn('[%s] "empty" in _connect is depricated and will be '
                        'ignored.' % self.get_step_name())
            set_out_connections.discard('empty')

        # For each parent step ...
        for parent in self.get_dependencies():
            if not parent.get_runs():
                raise UAPError('The step "%s" produces no output.' %
                        parent.get_step_name())
            logger.debug('Connecting "%s" to "%s".' %
                    (parent.get_step_name(), self.get_step_name()))
            # ... look for connection to add
            used_conns = cc.connect(parent, self, self._options['_connect'])
            if not used_conns:
                # ... or add connections with the same name.
                logger.debug('Parent "%s" not connected to child "%s". '
                        'Hence connecting equally named connections.'%
                         (parent.get_step_name(), self.get_step_name()))
                used_conns = cc.connect(parent, self)
            if not used_conns:
                raise UAPError('No connections could be made between '
                        '"%s" and its dependency "%s".' %
                        (self.get_step_name(), parent.get_step_name()))
            used_out_connections = used_out_connections.union(used_conns)

        # Check if all required connections are sattisfied.
        required_connections = self.get_in_connections(with_optional=False)
        missing = required_connections - cc.existing_connections
        if missing:
            logger.warn('_connect: The required connection %s of step '
                '"%s" is not satisfied. To remove this warning pass '
                'optional=True to the add_connection method in the step '
                'constructor __init__ of "%s".' %
                (missing, self.get_step_type(), self.get_step_type()))
            logger.warn('[Deprecation] Unmet required connections may trigger '
                'an error in future version of the UAP.')

        # Check if all set out connections were recognized.
        unrecognized = set_out_connections - used_out_connections
        if len(unrecognized) > 0:
            raise UAPError('For the following connections into step "%s" '
                    'no parent run could be found: %s.' %
                    (self.get_step_name(), list(unrecognized)))

        return cc


    def get_input_run_info_for_connection(self, in_key):
        if in_key[0:3] != 'in/':
            raise UAPError("in_key does not start with 'in/': %s" % in_key)
        if in_key not in self._connections:
            raise UAPError("Undeclared connection %s." % in_key)
        out_key = in_key.replace('in/', 'out/')
        out_keys = None
        allowed_steps = None
        if '_connect' in self._options:
            if in_key in self._options['_connect']:
                declaration = self._options['_connect'][in_key]
                if isinstance(declaration, str):
                    if '/' in declaration:
                        parts = declaration.split('/')
                        allowed_steps = set()
                        allowed_steps.add(parts[0])
                        out_key = 'out/' + parts[1]
                    else:
                        out_key = 'out/' + declaration
                elif isinstance(declaration, list):
                    for dec in declaration:
                        if isinstance(declaration, str):
                            if '/' in declaration:
                                parts = declaration.split('/')
                                if allowed_steps == None:
                                    allowed_steps = set()
                                allowed_steps.add(parts[0])
                                if out_keys == None:
                                    out_keys = list()
                                out_keys.append('out/' + parts[1])
                            else:
                                if out_keys == None:
                                    out_keys = list()
                                out_keys.append('out/' + declaration)
                else:
                    raise UAPError(
                        "Invalid _connect value: %s" % yaml.dump(declaration))

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
                result['counts'][key2] = (min if mkey == 'min' else max)\
                                         (result['counts'][key2], value)

        result['runs'] = dict()
        for step_name, step_info in self.get_input_runs().items():
            if allowed_steps is not None:
                if not step_name in allowed_steps:
                    continue
            for key in self.get_pipeline().get_step(step_name)._connections:
                if out_key == 'out/*' or out_key == key:
                    result['counts']['total_steps'] += 1
                    for run_id, run_info in step_info.items():
                        result['counts']['total_runs'] += 1
                        paths = run_info.get_output_files_abspath()[key.replace(
                            'out/', '')].keys()
                        result['counts']['total_files'] += len(paths)
                        if not run_id in result['runs']:
                            result['runs'][run_id] = dict()
                        result['runs'][run_id][step_name] = paths

                        steps_per_run = len(result['runs'][run_id])
                        update_min_max('steps_per_run', steps_per_run)

                        files_per_step_and_run = len(result['runs'][run_id]\
                                                     [step_name])
                        update_min_max('files_per_step_and_run',
                                       files_per_step_and_run)

                        files_per_run = 0
                        for _ in result['runs'][run_id].values():
                            files_per_run += len(_)
                        update_min_max('files_per_run', files_per_run)

        # check constraints, if any
        if in_key in self._connection_restrictions:
            for k, v in self._connection_restrictions[in_key].items():
                if result['counts'][k] != v:
                    raise UAPError("Connection constraint failed: %s/%s"
                                 "/%s should be %d but is %s."
                                 % (self, in_key, k, v,
                                    str(result['counts'][k])))

        return result


class AbstractSourceStep(AbstractStep):
    """
    A subclass all source steps inherit from and which distinguishes source
    steps from all real processing steps because they do not yield any tasks,
    because their "output files" are in fact files which are already there.

    Note that the name might be a bit misleading because this class only
    applies to source steps which 'serve' existing files. A step which has
    no input but produces input data for other steps and actually has to do
    something for it, on the other hand, would be a normal AbstractStep
    subclass because it produces tasks.
    """

    def __init__(self, pipeline):
        super(AbstractSourceStep, self).__init__(pipeline)
