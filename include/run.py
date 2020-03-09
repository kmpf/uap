import sys
from datetime import datetime, timedelta
import glob
import json
import fscache
from logging import getLogger
import os
import pwd
import random
import stat
import string
import tempfile
import platform
from deepdiff import DeepHash, DeepDiff
from collections import OrderedDict
import inspect
from functools import wraps

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
        self.fsc = fscache.FSCache()
        '''
        A cache.
        '''
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

    def cache(func):
        '''
        A decorator to cache a functions return value with self.fsc.
        '''
        function_name = [func.func_name, inspect.getargspec(func)]
        @wraps(func)
        def inner(self, *args, **kwargs):
            key = str(function_name + [args, kwargs])
            cache = self.fsc.cache
            if key in cache:
                result = cache[key]
            else:
                result = func(self, *args, **kwargs)
                cache[key] = result
            return result
        return inner


    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def reset_fsc(self):
        self.fsc.clear()

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

    def get_known_paths(self):
        return self._known_paths

    def add_known_paths(self, known_paths_dict):
        self._known_paths.update(known_paths_dict)

    def get_temp_paths(self):
        '''
        Returns a set of all temporary paths which belong to this run.
        '''
        return self._temp_paths

    def get_output_directory_du_jour_placeholder(self):
        '''
        Used to return a placeholder for the temporary output directory, which
        needed to be replaced by the actual temp directory inside the
        abstract_step.execute() method.
        '''
        raise UAPError("Using run.get_output_directory_du_jour_placeholder() "
                       "is deprecated. Just use the string '.' instead.")

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

    @cache
    def get_run_structure(self, commands=True):
        '''
        Returns a dictionary with all information known at
        run declaratuon time, relevant for its result
        and nothing more.

        Included are:
         - tool versions
         - commands and structure
         - output connections and files
         - parent run names and hashsum of their run_structure

        Should not include:
         - any absolute paths
         - parent hashes of source steps (absolute paths)
        '''

        cmd_by_eg = dict()
        cmd_by_eg['run_id'] = self.get_run_id()
        p = self.get_step().get_pipeline()

        # get tool version texts and paths
        tools = sorted(self.get_step()._tools.keys())
        cmd_by_eg['tool_versions'] = dict()
        tool_conf = p.config['tools']
        for tool in tools:
            if not tool_conf[tool]['ignore_version'] \
            and not p.args.no_tool_checks:
                tool_info = p.tool_versions[tool]
                real_tool_path = tool_info['used_path']
                response = tool_info['response'].replace(real_tool_path, tool)
                cmd_by_eg['tool_versions'][tool] = response

        # get output files
        cmd_by_eg['output'] = dict()
        for connection, files in sorted(self.get_output_files().items()):
            if files:
                cmd_by_eg['output'][connection] = files.keys()

        # get parent hash
        cmd_by_eg['parent hashes'] = dict()
        parents = self.get_parent_runs()
        for prun in parents:
            if isinstance(prun.get_step(), abst.AbstractSourceStep):
                continue
            task_id = '%s/%s' % (prun.get_step().get_step_name(), prun.get_run_id())
            hashsum = misc.str_to_sha256(json.dumps(prun.get_run_structure(),
                    sort_keys=True))
            cmd_by_eg['parent hashes'][task_id] = hashsum

        if not commands:
            return cmd_by_eg

        # get commands
        for eg_count, exec_group in enumerate(self.get_exec_groups()):
            eg_name = 'execution group %d' % eg_count
            cmd_by_eg[eg_name] = dict()
            procs = exec_group.get_pipes_and_commands(sort=True)
            for pipe_count, poc in enumerate(procs):
                # for each pipe or command (poc)
                # check if it is a pipeline ...
                if isinstance(poc, pipeline_info.PipelineInfo):
                    cmd_by_eg[eg_name]['pipe %s' % pipe_count] = \
                            poc.get_command_string(replace_path=True)
                # ... or a command
                elif isinstance(poc, command_info.CommandInfo):
                    cmd_by_eg[eg_name]['command %s' % pipe_count] = \
                            poc.get_command_string(replace_path=True)

        return cmd_by_eg

    def get_changes(self):
        anno_data = self.written_anno_data()
        if not anno_data:
            return {'Error':'Missing annotation file %s' %
                    self.get_annotation_path()}
        old_struct = anno_data['run']['structure']
        new_struct = self.get_run_structure()
        return DeepDiff(old_struct, new_struct)

    def dependencies(self):
        """
        Returns a dict with a set of input files for each output file.
        """
        connections = self.get_output_files_abspath()
        dependencies = dict()
        for deps in connections.values():
            for out_file, input_files in deps.items():
                if not input_files or input_files == [None]:
                    continue
                dependencies.setdefault(out_file, set())
                dependencies[out_file].update(set(input_files))
        return dependencies

    def file_changes(self, do_hash=False, report_correct=False):
        anno_data = self.written_anno_data()
        if not anno_data:
            raise StopIteration
        new_dest = self.get_step().get_pipeline().config['destination_path']
        old_dest = anno_data['config']['destination_path']
        end_time = anno_data['end_time']
        p = self.get_step().get_pipeline()
        is_volatile = self.get_step().is_volatile()
        for path, input_files in self.dependencies().items():

            # is it logged in the annotation file
            old_path = path.replace(new_dest, old_dest)
            if old_path not in anno_data['run']['known_paths']:
                yield '%s not logged in the annotation file %s' % \
                      (old_path, os.path.basename(self.get_annotation_path()))
                continue
            meta_data = anno_data['run']['known_paths'][old_path]

            # existence and volatility
            v_path = path + abst.AbstractStep.VOLATILE_SUFFIX
            if not self.fsc.exists(path):
                if not is_volatile or not self.fsc.exists(v_path):
                    yield '%s is missing' % path
                    continue
                elif is_volatile:
                    path = v_path
                else:
                    # is marked volatile but not volatilized yet
                    is_volatile = False

            # modification time
            change_str = ''
            new_mtime = datetime.fromtimestamp(self.fsc.getmtime(path))
            old_mtime = meta_data['modification time']
            diff = new_mtime - old_mtime
            if diff > timedelta(seconds=1):
                change_str = ' and modification date after %s' % old_mtime

            # chaged dependencies
            has_changed_deps = False
            for in_file in input_files:
                parent_task = p.get_task_for_file(in_file)
                if parent_task:
                    parent_volitile = parent_task.step.is_volatile()
                    parent_fsc = parent_task.get_run().fsc
                else:
                    # the parent was probably a source step
                    parent_volitile = False
                    parent_fsc = self.fsc
                v_in = in_file + abst.AbstractStep.VOLATILE_SUFFIX
                if not parent_fsc.exists(in_file):
                    if not parent_volitile or not parent_fsc.exists(v_in):
                        has_changed_deps = True
                        yield 'input file %s is missing' % in_file
                        continue
                    else:
                        in_file = v_in
                if parent_fsc.getmtime(in_file) > \
                self.fsc.getmtime(path):
                    has_changed_deps = True
                    yield 'input file %s was modified' % in_file
            if has_changed_deps:
                if change_str:
                    change_str = ', has changed input' + change_str
                else:
                    change_str = ' and has changed input'

            # stop here if volatile
            if is_volatile and change_str:
                yield path + ' is volatilized' + change_str
                continue
            elif is_volatile:
                continue

            # size changes
            old_size = meta_data['size']
            new_size = self.fsc.getsize(path)
            if new_size != old_size:
                yield '%s size changed from %s B to %s B%s' % \
                      (path, old_size, new_size, change_str)
                continue

            # hash sum
            if do_hash is True:
                old_hash = meta_data['sha256']
                new_hash = self.fsc.sha256sum_of(path)
                if new_hash != old_hash:
                    yield '%s sha256sum changed from %s to %s%s' % \
                          (path, old_hash, new_hash, change_str)
                    continue
                elif change_str or report_correct is True:
                    yield '%s sha256sum is correct%s' % \
                          (path, change_str)
                    continue

            # what changed
            elif change_str:
                if change_str.startswith(' and'):
                    change_str = change_str[len(' and'):]
                elif change_str.startswith(' ,'):
                    change_str = change_str[len(' ,'):]
                yield path + change_str

    @cache
    def get_state(self, do_hash=False, reset=False):

        states = self.get_step().get_pipeline().states
        if isinstance(self.get_step(), abst.AbstractSourceStep):
            return states.FINISHED
        ex_ping_file = self.get_executing_ping_file()
        if self.fsc.exists(ex_ping_file):
            logger.debug('Found execution ping file: %s' % ex_ping_file)
            if self.is_stale():
                return states.BAD
            return states.EXECUTING
        qu_ping_file = self.get_queued_ping_file()
        if self.fsc.exists(qu_ping_file):
            logger.debug('Found queue ping file: %s' % qu_ping_file)
            return states.QUEUED
        if self.fsc.exists(self.get_queued_ping_file() + '.bad'):
            return states.BAD

        anno_data = self.written_anno_data()
        if anno_data:
            if anno_data.get('run', dict()).get('error'):
                return states.BAD
            if self.get_changes():
                return states.CHANGED

        has_volitile_parent = False
        for parent in self.get_parent_runs():
            pstate = parent.get_state(do_hash=do_hash)
            if pstate not in [states.FINISHED, states.VOLATILIZED]:
                return states.WAITING
            if pstate == states.VOLATILIZED:
                has_volitile_parent = True

        output_files = [(out_file, input_files)
                for files in self.get_output_files_abspath().values()
                for out_file, input_files in files.items()]
        all_exist = all(self.fsc.exists(out_file)
                        for out_file, _ in output_files)
        is_volatilized = False
        if not all_exist and self.get_step().is_volatile():
            all_exist = all(self.fsc.exists(out_file
                            + abst.AbstractStep.VOLATILE_SUFFIX)
                            for out_file, _ in output_files)
            is_volatilized = True

        if not all_exist:
            if has_volitile_parent:
                return states.WAITING
            return states.READY
        else:
            if not anno_data:
                return states.CHANGED
            for bad_file in self.file_changes(do_hash=do_hash):
                if bad_file:
                    return states.CHANGED
            if is_volatilized:
                return states.VOLATILIZED
            return states.FINISHED

    def get_parent_runs(self):
        """
        Returns the parent runs.
        """
        p = self.get_step().get_pipeline()
        task_id = '%s/%s' % (self.get_step(), self.get_run_id())
        input_files = set()
        if task_id in p.input_files_for_task_id:
            input_files = p.input_files_for_task_id[task_id]
        parents = set()
        # Only source steps do have empty strings in the input files list
        # so we can safely exclude them here
        for inpath in [x for x in input_files if x != '']:
            task_id = p.task_id_for_output_file[inpath]
            if task_id in p.task_for_task_id:
                parents.add(p.task_for_task_id[task_id].get_run())
        return parents

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
        else:
            out_path = os.path.normpath(out_path)
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

    def add_temporary_file(self, prefix = 'temp', suffix = '', designation = None):
        '''
        Returns the name of a temporary file.
        '''
        count = len(self._temp_paths)
        count = 0
        while True:
            temp_name = prefix + '-' + str(count) + suffix
            if not temp_name in self._temp_paths:
                break
            else:
                count += 1

        logger.debug("Temporary file (#%s) in run %s: %s" %
              (len(self._temp_paths) + 1, self.get_run_id(), temp_name) )

        # _known_paths dict is logged
        known_paths = dict()
        known_paths[temp_name] = {
            'label': os.path.basename(temp_name),
            'designation': designation,
            'type': ''
        }
        self.add_known_paths(known_paths)
        # _temp_paths set contains all temporary files which are going to be
        # deleted
        self._temp_paths.add(temp_name)
        return temp_name

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
        for tmp_file in self._temp_paths:
            # Check file type
            if not os.path.exists(tmp_file):
                logger.debug("Set %s 'type' info to 'missing'" % tmp_file)
                self._known_paths[tmp_file]['type'] = 'missing'
                continue
            pathmode = os.stat(tmp_file).st_mode
            isdir = stat.S_ISDIR(pathmode) != 0
            isfile = stat.S_ISREG(pathmode) != 0
            isfifo = stat.S_ISFIFO(pathmode) != 0
            # Update 'type' value
            if tmp_file in self._known_paths.keys():
                if isfile:
                    logger.debug("Set %s 'type' info to 'file'" % tmp_file)
                    self._known_paths[tmp_file]['type'] = 'file'
                elif isdir:
                    logger.debug("Set %s 'type' info to 'directory'" % tmp_file)
                    self._known_paths[tmp_file]['type'] = 'directory'
                elif isfifo:
                    logger.debug("Set %s 'type' info to 'fifo'" % tmp_file)
                    self._known_paths[tmp_file]['type'] = 'fifo'
            if os.path.isdir(tmp_file) and isdir:
                try:
                    logger.info("Now deleting directory: %s" % tmp_file)
                    os.rmdir(tmp_file)
                except OSError as e:
                    logger.error("errno: %s" % e.errno)
                    logger.error("strerror: %s" % e.strerror)
                    logger.error("filename: %s" % e.filename)
                    pass
            else:
                try:
                    logger.info("Now deleting: %s" % tmp_file)
                    os.unlink(tmp_file)
                except OSError as e:
                    pass

    def add_empty_output_connection(self, tag):
        '''
        An empty output connection has 'None' as output file and 'None' as input
        file.
        '''
        logger.warn('[Deprecation] %s: add_empty_output_connection is depricated. '
                'Please make the connection "out/%s" optional and do not add '
                'anything instead.' % (self.get_step().get_step_type(), tag))
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
        raise UAPError("Sorry, your output '%s' file couldn't be found in"
                     "the dictionary: %s." % (out_path, temp))

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
                directory = self.get_output_directory()
                full_path = out_path
                try:
                    if directory and out_path and directory != '.':
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
        result = OrderedDict([
            ('step_type', self.get_step().get_step_type()),
            ('step_name', self.get_step().get_step_name()),
            ('run_id', self._run_id),
            ('output_directory', self.get_output_directory()),
            ('annotation_file', self.get_annotation_path()),
            ('output_files', self._output_files),
            ('private_info', self._private_info),
            ('public_info', self._public_info)
        ])
        result.update(self.get_run_structure(commands=False))
        del result['tool_versions']
        del result['output']
        anno = self.written_anno_data()
        if anno and anno.get('run'):
            result['run'] = OrderedDict([
                ('start_time', anno['start_time']),
                ('end_time', anno['end_time'])
            ])
            for key, value in anno['run'].items():
                if key in result.keys() + ['known_paths', 'structure']:
                    continue
                result['run'][key] = value
        else:
            result['run'] = 'not run yet'
        return result

    @cache
    def written_anno_data(self):
        anno_file = self.get_annotation_path()
        try:
            with open(anno_file, 'r') as fl:
                return yaml.load(fl, Loader=yaml.FullLoader)
        except IOError as e:
            if not os.path.exists(anno_file):
                return False
            else:
                logger.warn('The annotation file "%s" could not be read.'
                            % anno_file)
        return None

    def write_annotation_file(self, path=None, error=None, job_id=None):
        '''
        Write the YAML annotation after a successful or failed run. The
        annotation can later be used to render the process graph.
        '''
        if path is None:
            path = self.get_output_directory()

        # now write the annotation
        log = {}
        log['pid'] = os.getpid()
        log['step'] = {}
        log['step']['options'] = self.get_step().get_options()
        log['step']['name'] = self.get_step().get_step_name()
        log['step']['type'] = self.get_step().get_step_type()
        # if a submit script was used ...
        script = self.get_step().get_submit_script_file()
        script_path = os.path.relpath(script, path)
        if os.path.exists(script):
            log['step']['submit_script'] = script_path
        log['step']['cores'] = self.get_step().get_cores()
        log['run'] = {}
        log['run']['run_id'] = self.get_run_id()
        log['run']['output_directory'] = self.get_output_directory()
        log['run']['private_info'] = self._private_info
        log['run']['public_info'] = self._public_info
        log['run']['temp_directory'] = self.get_temp_output_directory()
        # if a run submit script was used ...
        if os.path.exists(self.get_submit_script_file()):
            # ... read it and store it ...
            with open(self.get_submit_script_file(), 'r') as f:
                log['run']['submit_script'] = f.read()
            # ... finally delete it
            os.unlink(self.get_submit_script_file())
        log['run']['known_paths'] = self.get_known_paths()
        log['run']['structure'] = self.get_run_structure()
        log['run']['hostname'] = platform.node()
        log['run']['platform'] = platform.platform()
        log['run']['user'] = pwd.getpwuid(os.getuid())[0]
        if error is not None:
            log['run']['error'] = error
        if job_id:
            log['run']['cluster job id'] = job_id
        else:
            try:
                with open(self.get_queued_ping_file(), 'r') as buff:
                    info = yaml.load(buff, Loader=yaml.FullLoader)
                log['run']['cluster job id'] = info['job_id']
            except (IOError, KeyError):
                pass

        p = self.get_step().get_pipeline()
        log['config'] = p.config
        if p.args.no_tool_checks:
            logger.warn('Writing annotation file without tool checks.')
            log['tool_versions'] = 'deactivated with --no-tool-checks'
        else:
            log['tool_versions'] = {}
            for tool in self.get_step()._tools.keys():
                log['tool_versions'][tool] = p.tool_versions[tool]
        log['pipeline_log'] = self.get_step()._pipeline_log
        log['start_time'] = self.get_step().start_time
        log['end_time'] = self.get_step().end_time



        log['uap_version'] = p.args.uap_version
        log['git_tag'] = p.git_tag
        log['git_diff'] = p.git_diff
        log['git_version'] = p.git_version


        if p.caught_signal is not None:
            log['signal'] = p.caught_signal

        annotation_yaml = yaml.dump(log, default_flow_style = False)
        annotation_path = self.get_annotation_path(path)

        # overwrite the annotation if it already exists
        with open(annotation_path, 'w') as f:
            f.write(annotation_yaml)

        return annotation_path

    def get_annotation_path(self, path=None):
        if path is None:
            path = self.get_output_directory()
        annotation_path = os.path.join(
            path, "%s-annotation.yaml" % self.get_run_id()
        )
        return annotation_path

    def is_stale(self, exec_ping_file=None):
        """
        Returns time of inactivity if the ping file exists and is stale.
        """
        if exec_ping_file is None:
            exec_ping_file = self.get_executing_ping_file()
        if self.fsc.exists(exec_ping_file):
            last_activity = datetime.fromtimestamp(
                    self.fsc.getmtime(exec_ping_file))
            now = datetime.now()
            inactivity = now - last_activity
            if inactivity.total_seconds() > \
                    abst.AbstractStep.PING_TIMEOUT:
                        return inactivity
        return False
