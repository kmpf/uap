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
import platform
import subprocess

import yaml

import abstract_step as abst
import command as command_info
import exec_group
import pipeline_info
import misc
from uaperrors import UAPError

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
            raise UAPError("Error: A run ID must not contain a slash: %s." % run_id)
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
        self._input_files = set()
        self._output_files = dict()
        for out_connection in self._step.get_out_connections(with_optional=False):
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
        self._temp_paths = set()
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
            raise KeyError("Connection %s not declared for step %s" %
                         (connection, self.get_step()))

    @property
    def annotation_file(self):
        afl = glob.glob( os.path.join(
            self.get_output_directory(),
            ".%s-annotation-*.yaml" % self.get_run_id()))
        if not afl:
            print(self.get_output_directory())
            return ""
        elif len(afl) != 1:
            raise StandardError("Found multiple annotation files: %s" %
                                ", ".join(afl))
        elif os.path.isfile(afl[0]):
            return afl[0]

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
        if self._submit_script == None:
            self._submit_script = os.path.join(
                self.get_output_directory(),
                ".submit-%s-%s.sh" % (self.get_step().get_step_name(),
                                      self.get_run_id())
            )
        return self._submit_script

    def is_source(self):
        return True if isinstance(self._step, abst.AbstractSourceStep) else False

    def replace_output_dir_du_jour(func):
        def inner(self, *args, **kwargs):
            # Collect info to replace du_jour placeholder with temp_out_dir
            step = self.get_step()
            placeholder = self.get_output_directory_du_jour_placeholder()
            temp_out_dir = self.get_output_directory_du_jour()

            value = None
            ret_value = func(self, *args, **kwargs)
            # If currently calling AbstractStep.runs() do nothing
            if ret_value is None:
                return(None)
            elif temp_out_dir is None:
                return(ret_value)
            elif isinstance(ret_value, list) or isinstance(ret_value, set):
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
                raise UAPError("Function %s does not return list or string object"
                             % func.__class__.__name__)
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
        Returns a set of all temporary paths which belong to this run.
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
           this will return the current directory
         - otherwise, it will return the real output directory
        '''
        if self.get_step()._state == abst.AbstractStep.states.DECLARING:
            return None
        elif self.get_step()._state == abst.AbstractStep.states.EXECUTING:
            return '.'
        else:
            return self.get_output_directory()

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

    def get_run_structure(self):
        '''
        Creates a dictionary with the structure of commands to
        be executed and the used tool versions.

        This causes runs to be marked for rerunning if the commands to be
        executed change.
        '''

        step = self.get_step()
        # Store step state
        previous_state = step._state
        # Set step state to DECLARING to avoid circular dependencies
        step._state = abst.AbstractStep.states.DECLARING

        cmd_by_eg = dict()
        # get tool version texts
        tools = sorted(step._tools.keys())
        cmd_by_eg['tool_versions'] = dict()
        tool_paths = dict()
        tool_conf = step.get_pipeline().config['tools']
        for tool in tools:
            tool_info = step.get_pipeline().tool_versions[tool]
            if tool != tool_conf[tool]['path']:
                tool_paths[tool] = tool_conf[tool]['path']
            if tool_conf[tool]['ignore_version'] is not True:
                real_tool_path = tool_info['used_path']
                response = tool_info['response'].replace(real_tool_path, tool)
                cmd_by_eg['tool_versions'][tool] = response

        # get commands
        for eg_count, exec_group in enumerate(self.get_exec_groups()):
            eg_name = 'execution group %d' % eg_count
            cmd_by_eg[eg_name] = dict()
            for pipe_count, poc in enumerate(exec_group.get_pipes_and_commands()):
                # for each pipe or command (poc)
                # check if it is a pipeline ...
                if isinstance(poc, pipeline_info.PipelineInfo):
                    pipe_count += 1
                    cmd_by_eg[eg_name]['pipe %s' % pipe_count] = list()
                    for command in poc.get_commands():
                        cmd_list = command.get_command()
                        # replace tool paths by tool names
                        for tool, path in tool_paths.items():
                            cmd_list = [element.replace(path, tool)
                                    for element in cmd_list]
                        cmd_by_eg[eg_name]['pipe %s' % pipe_count].append(
                                subprocess.list2cmdline(cmd_list))
                # ... or a command
                elif isinstance(poc, command_info.CommandInfo):
                    cmd_by_eg[eg_name]['command %s' % pipe_count] = \
                            subprocess.list2cmdline(poc.get_command())

        # Set step state back to original state
        step._state = previous_state
        return cmd_by_eg

    def get_output_directory(self):
        '''
        Returns the final output directory.
        '''
        return os.path.join(
            self.get_step().get_output_directory(),
            self.get_run_id()
        )

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
            raise UAPError(
                "You're trying to overwrite private info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._private_info[key]))
        self._private_info[key] = value

    def add_public_info(self, key, value):
        '''
        Add public information to a run. For example, a FASTQ reader may store
        the index barcode here for subsequent steps to query via
        ``AbstractStep.find_upstream_info()``.
        '''
        if key in self._public_info and value != self._public_info[key]:
            raise UAPError(
                "You're trying to overwrite public info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._public_info[key]))
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
            raise UAPError("The key %s doesn't exist yet as public info."
                         "Please use add_public_info(%s, %s)"
                         % (key, key, value))
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
            raise UAPError("The declared output file path contains "
                         "directory separator: %s." % out_path)
        elif isinstance(self._step, abst.AbstractSourceStep):
            out_path = os.path.abspath(out_path)
        # make sure tag was declared with an outgoing connection
        if 'out/' + tag not in self._output_files:
            self.add_out_connection('out/' + tag)

        out_connection = self.get_out_connection(tag)

        if out_path in self.get_output_files_for_out_connection(out_connection):
            raise UAPError(
                "You're trying to re-add an output file which has already "
                "been declared: %s." % out_path)

        if not isinstance(in_paths, list):
            raise UAPError("Input paths (%s) is not a list." % in_paths)

        if None in in_paths:
            raise UAPError(
                "There is a NoneType element in input paths (%s) for output "
                "file (%s)" % (in_paths, out_path))

        if out_path == None:
            raise UAPError(
                "Trying to add NoneType element as output file for input paths "
                ": %s" % in_paths)

        self._input_files.union(set(in_paths))
        logger.debug('Adding files %s as for connection %s in %s for run %s.' %
                (out_path, out_connection, str(self.get_step()), self.get_run_id()))
        self._output_files[out_connection][out_path] = in_paths
        return out_path

    @replace_output_dir_du_jour
    def add_temporary_file(self, prefix = 'temp', suffix = '', designation = None):
        '''
        Returns the name of a temporary file (created by tempfile library).
        Name and output directory placeholder are concatenated. The concatenated
        string is returned and stored in a list. The placeholder is immediately
        properly adjusted by @replace_output_dir_du_jour.
        '''
        count = len(self._temp_paths)
        temp_placeholder = str()

        count = 0
        while True:
            temp_name = prefix + '-' + str(count) + suffix
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
        # _temp_paths set contains all temporary files which are going to be
        # deleted
        self._temp_paths.add(temp_placeholder)
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
        possible. The list elements are removed in LIFO order.
        Also, self._known_paths 'type' info is updated here.
        NOTE: Included additional stat checks to detect FIFOs as well as other
        special files.
        '''
        for _ in self.get_temp_paths()[::-1]:
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
        logger.warn('[Deprecation] add_empty_output_connection is depricated. '
                'Please make the connection "out/%s" optional and do not add '
                'anything instead.' % tag)
        # make sure tag was declared with an outgoing connection
        if 'out/' + tag not in self._step.get_out_connections():
            raise UAPError("Invalid output_file tag '%s' in %s. "
                         "You might want to add self.add_connection('out/%s') "
                         "to the constructor of %s."
                         % (tag, str(self._step), tag, self._step.__module__))
        try:
            out_connection = self.get_out_connection(tag)
        except KeyError:
            out_connection = self.add_out_connection(tag)

        if None in self._output_files[out_connection]:
            raise UAPError(
                "You're trying to re-declare %s as an empty output connection "
                % out_connection)

        self._output_files[out_connection][None] = None

    def add_out_connection(self, out_connection):
        if not out_connection.startswith('out/'):
            out_connection = 'out/' + out_connection
        if out_connection not in self._step.get_out_connections():
            raise UAPError("Invalid output connection '%s' in %s. "
                         "You might want to add self.add_connection('%s') "
                         "to the constructor of %s."
                         % (out_connection, str(self._step), out_connection,
                                 self._step.__module__))
        logger.debug('Adding %s to %s in run %s.' %
                (out_connection, str(self.get_step()), self.get_run_id()))
        self._output_files[out_connection] = dict()
        return out_connection

    def get_input_files_for_output_file(self, output_file):
        for connection in self.get_out_connections():
            if output_file in \
                self.get_output_files_for_out_connection(connection):
                    return self._output_files[connection][output_file]

    def get_input_files_for_output_file_abspath(self, output_file):
        for connection in self.get_out_connections():
            if output_file in \
               self.get_output_files_abspath_for_out_connection(connection):
                return self.get_output_files_abspath()[connection]\
                    [output_file]

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
            raise UAPError("More than one output file declared for out/%s."
                         % annotation)
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
        raise UAPError("Sorry, your output '%s' file couldn't be found in"
                     "the dictionary: %s." % (out_path, temp))

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
        # if a submit script was used ...
        script = self.get_step().get_submit_script_file()
        if os.path.exists(script):
            # ... read it and store it ...
            with open(script, 'r') as f:
                log['step']['submit_script'] = f.read()
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
        log['run']['structure'] = self.get_run_structure()
        log['config'] = self.get_step().get_pipeline().config

        log['tool_versions'] = {}
        for tool in self.get_step()._tools.keys():
            log['tool_versions'][tool] = self.get_step().get_pipeline()\
                                                        .tool_versions[tool]
        log['pipeline_log'] = self.get_step()._pipeline_log
        log['start_time'] = self.get_step().start_time
        log['end_time'] = self.get_step().end_time



        log['git_tag'] = self.get_step().get_pipeline().git_tag
        log['git_diff'] = self.get_step().get_pipeline().git_diff
        log['git_version'] = self.get_step().get_pipeline().git_version
        log['system'] = {}
        log['system']['hostname'] = platform.node()
        log['system']['platform'] = platform.platform()


        if self.get_step().get_pipeline().caught_signal is not None:
            log['signal'] = self.get_step().get_pipeline().caught_signal

        annotation_yaml = yaml.dump(log, default_flow_style = False)
        annotation_path = os.path.join(
            path, "%s-annotation.yaml" % self.get_run_id()
        )

        # overwrite the annotation if it already exists
        with open(annotation_path, 'w') as f:
            f.write(annotation_yaml)

        return annotation_path, annotation_yaml
