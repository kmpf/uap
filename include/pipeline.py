import base64
import copy
import datetime
import json
import signal
from logging import getLogger
from operator import itemgetter
import os
import re
import subprocess
import sys
import yaml
import multiprocessing
import traceback

import abstract_step
import misc
import task as task_module
from uaperrors import UAPError


logger = getLogger("uap_logger")


class ConfigurationException(Exception):
    """an exception class for reporting configuration errors"""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def exec_pre_post_calls(tool_id, info_key, info_command,
                        tool_check_info):
    if isinstance(info_command, str):
        info_command = [info_command]
    for command in info_command:
        if isinstance(command, str):
            command = command.split()
        for argument in command:
            if not isinstance(argument, str):
                raise UAPError(
                    "The command to be launched '%s' contains non-string "
                    "argument '%s'. Therefore the command will fail. Please "
                    "fix this type issue." % (command, argument))
        logger.info("Executing command: %s" % " ".join(command))
        try:
            proc = subprocess.Popen(
                command,
                stdin = None,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                close_fds = True)

        except OSError as e:
            raise UAPError(
                "Error while executing '%s' for %s: %s "
                "Error no.: %s Error message: %s" %
                (info_key, tool_id,
                 " ".join(command), e.errno, e.strerror)
            )

        command_call = info_key
        command_exit_code = '%s-exit-code' % info_key
        command_response = '%s-respone' % info_key
        (output, error) = proc.communicate()
        if info_key in ['module_load', 'module_unload']:
            logger.info("Try '%s' for '%s': %s" % (
                info_key, tool_id, " ".join(command))
            )
            try:
                exec output
            except NameError as e:
                msg = "Error while loading module '%s': \n%s"
                raise UAPError(msg % (tool_id, error))

            tool_check_info.update({
                command_call : (' '.join(command)).strip(),
                command_exit_code : proc.returncode
            })
            sys.stderr.write(error)
            sys.stderr.flush()
        else:
            tool_check_info.update({
                command_call : (' '.join(command)).strip(),
                command_exit_code : proc.returncode,
                command_response : (output + error)
            })

    return tool_check_info

def check_tool(args):
    '''
    A top-level function to be used a multiprocessing pool to retrieve tool
    information in parallel.
    '''
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        tool_id, info = args
        tool_check_info = dict()

        # Load module(s) and execute command if configured
        for pre_cmd in (x for x in ('module_load', 'pre_command')
                         if x in info):
            tool_check_info = exec_pre_post_calls(
                tool_id, pre_cmd, info[pre_cmd], tool_check_info)

        # Execute command to check if tool is available
        command = info['path']
        if isinstance(command, str):
            command = [command]
        elif not isinstance(command, list):
            raise TypeError('Unsupported format for path of tool "%s": %s' %
                    (tool_id, info['path']))
        command.append(info['get_version'])
        logger.info("Executing command: %s" % " ".join(command))
        try:
            proc = subprocess.Popen(
                command,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                close_fds = True)
            proc.stdin.close()
        except OSError as e:
            raise UAPError("Error while checking Tool %s "
                         "Error no.: %s Error message: %s\ncommand: %s "
                         "\nSTDOUT-ERR: %s\n" %
                         (info['path'], e.errno, e.strerror, command,
                                 subprocess.PIPE))
        proc.wait()
        exit_code = None
        exit_code = proc.returncode
        tool_check_info.update({
            'command': (' '.join(command)).strip(),
            'exit_code': exit_code,
            'response': (proc.stdout.read() + proc.stderr.read()).strip()
        })
        expected_exit_code = info['exit_code']
        if exit_code != expected_exit_code:
            raise UAPError(
                "%s: Tool check failed for %s: %s - exit code is: %d "
                "(expected %d) (response %s)"
                % (self.config_name, tool_id, ' '.join(command),
                   exit_code, expected_exit_code, tool_check_info['response'])
            )
        # Execute clean-up command (if configured)
        for info_key in (x for x in ('module_unload', 'post_command')
                         if x in info):
            tool_check_info = exec_pre_post_calls(
                tool_id, info_key, info[info_key], tool_check_info)
    except Exception as e:
        logger.error(traceback.format_exc())
        raise
    return tool_id, tool_check_info

class Pipeline(object):
    '''
    The Pipeline class represents the entire processing pipeline which is defined
    and configured via the configuration file config.yaml.

    Individual steps may be defined in a tree, and their combination with samples
    as generated by one or more source leads to an array of tasks.
    '''

    states = misc.Enum(['WAITING', 'READY', 'QUEUED', 'EXECUTING', 'FINISHED'])
    '''
    Possible states a task can be in.
    '''

    def __init__(self, **kwargs):
        self.caught_signal = None
        self._cluster_type = None
        self.git_version = None
        self.git_diff = None
        self.git_tag = None

        '''use git diff to determine any changes in git
        directory if git is available
        '''
        command = ['git', '--version']
        try:

            self.git_version = subprocess.check_output(command).strip()

        except subprocess.CalledProcessError:
            logger.warn("""Execution of %s failed. Git seems to be
                         unavailable. Continue anyways""" % " ".join(command))

        if self.git_version:
            command = ['git', 'diff']
            try:
                self.git_diff = subprocess.check_output(command)
            except subprocess.CalledProcessError:
                logger.error("Execution of %s failed." % " ".join(command))

            command = ['git', 'describe', '--all', '--long']
            try:
                self.git_tag = subprocess.check_output(command).strip()
            except subprocess.CalledProcessError:
                logger.error("Execution of %s failed." % " ".join(command))

            if self.git_diff != '':
                logger.warn('THE GIT REPOSITORY HAS UNCOMMITED CHANGES!')


        """
        check if we got passed an 'arguments' parameter
        this parameter should contain a argparse.Namespace object
        """
        self.args = None
        if 'arguments' in kwargs:
            self.args = kwargs['arguments']


        '''
        Absolute path to the directory of the uap executable.
        It is used to circumvent path issues.
        '''
        self._uap_path = self.args.uap_path


        '''
        The cluster type to be used (must be one of the keys specified in
        cluster_config).
        '''
        self._cluster_config_path = os.path.join(
            self._uap_path, 'cluster/cluster-specific-commands.yaml')
        with open(self._cluster_config_path, 'r') as cluster_config_file:
            self._cluster_config = yaml.load(cluster_config_file, Loader=yaml.FullLoader)

        try:
            # set cluster type
            if self.args.cluster == 'auto':
                self.set_cluster_type(self.autodetect_cluster_type())
            else:
                self.set_cluster_type(self.args.cluster)
        except AttributeError:
            # cluster type is not an applicable parameter here, and that's fine
            # (we're probably in run-locally.py)
            pass

        self._start_working_dir = os.getcwd()
        '''
        User working directory.
        '''

        self.config_path, self.config_name = os.path.split(self.args.config.name)
        '''
        Name of the YAML configuration file
        '''

        self._config_path = os.path.abspath(self.config_path)
        '''
        Path of the YAML configuration file
        '''
        os.chdir(self._config_path)

        self.config = dict()
        '''
        Dictionary representation of configuration YAML file.
        '''

        self.steps = dict()
        '''
        This dict stores step objects by their name. Each step knows his
        dependencies.
        '''

        self.topological_step_order = list()
        '''
        List with topologically ordered steps.
        '''

        self.file_dependencies = dict()
        '''
        This dict stores file dependencies within this pipeline, but regardless
        of step, output file tag or run ID. This dict has, for all output
        files generated by the pipeline, a set of input files that output
        file depends on.
        '''

        self.file_dependencies_reverse = dict()
        '''

        This dict stores file dependencies within this pipeline, but regardless
        of step, output file tag or run ID. This dict has, for all input
        files required by the pipeline, a set of output files which are generated
        using this input file.
        '''

        self.task_id_for_output_file = dict()
        '''
        This dict stores a task ID for every output file created by the pipeline.
        '''

        self.task_ids_for_input_file = dict()
        '''
        This dict stores a set of task IDs for every input file used in the
        pipeline.
        '''

        self.input_files_for_task_id = dict()
        '''
        This dict stores a set of input files for every task id in the pipeline.
        '''

        self.output_files_for_task_id = dict()
        '''
        This dict stores a set of output files for every task id in the pipeline.
        '''

        self.task_for_task_id = dict()
        '''
        This dict stores task objects by task IDs.
        '''

        self.all_tasks_topologically_sorted = list()
        '''
        List of all tasks in topological order.
        '''

        self.tasks_in_step = dict()
        '''
        This dict stores tasks per step name.
        '''

        self.known_config_keys = set(['destination_path', 'constants', 'cluster',
                'steps', 'lmod', 'tools'])
        '''
        A set of accepted keys in the config.
        '''

        self.read_config(self.args.config)
        self.setup_lmod()

        self.tool_versions = {}
        self.check_tools()

        # collect all tasks
        for step_name in self.topological_step_order:
            step = self.get_step(step_name)
            self.tasks_in_step[step_name] = list()
            logger.debug("Collect now all tasks for step: %s" % step)
            for run_index, run_id in enumerate(misc.natsorted(step.get_run_ids())):
                task = task_module.Task(self, step, run_id, run_index)
                # if any run of a step contains an exec_groups,
                # the task (step/run) is added to the task list
                run = step.get_run(run_id)
                logger.debug("Step: %s, Run: %s" % (step, run_id))
                run_has_exec_groups = False
                if len(run.get_exec_groups()) > 0:
                    run_has_exec_groups = True
                if run_has_exec_groups:
                    logger.debug("Task: %s" % task)
                    self.all_tasks_topologically_sorted.append(task)
                    self.tasks_in_step[step_name].append(task)
                # Fail if multiple tasks with the same name exist
                if str(task) in self.task_for_task_id:
                    raise UAPError("%s: Duplicate task ID %s." %
                                 (self.config_name, str(task)))
                self.task_for_task_id[str(task)] = task

    def get_uap_path(self):
        return self._uap_path

    def get_config_path(self):
        '''
        Working directory of the pipeline that coinsides with the config path.
        '''
        return self._config_path

    def get_cluster_config(self):
        return self._cluster_config

    def get_steps(self):
        return self.steps

    def get_step(self, step_name):
        return self.steps[step_name]

    # read configuration and make sure it's good
    def read_config(self, config_file):
        # yaml.load works fine, even for duplicate dictionary keys (WTF)
        self.config = yaml.load(config_file, Loader=yaml.FullLoader)

        for key in self.config.keys():
            if key not in self.known_config_keys:
                raise UAPError('The key "%s" set in "%s" is unknown.' %
                        (key, self.args.config.name))

        if not 'id' in self.config:
            self.config['id'] = self.config_name

        if 'lmod' not in self.config or self.config['lmod'] is None:
            self.config['lmod'] = dict()
        if os.environ.has_key('LMOD_CMD'):
            self.config['lmod'].setdefault('path', os.environ['LMOD_CMD'])
        if os.environ.has_key('MODULEPATH'):
            self.config['lmod'].setdefault('module_path', os.environ['MODULEPATH'])
        for key in ('path', 'module_path'):
            if key not in self.config['lmod']:
                raise UAPError('lmod is not loaded or misses the key %s' % key)

        if not 'tools' in self.config or not isinstance(self.config['tools'], dict):
            self.config['tools'] = dict()
        for tool, args in self.config['tools'].items():
            if args is None or len(args)==0:
                self.config['tools'][tool] = dict()
            self.config['tools'][tool].setdefault('path', tool)
            self.config['tools'][tool].setdefault('get_version', '--version')
            self.config['tools'][tool].setdefault('exit_code', 0)
            self.config['tools'][tool].setdefault('ignore_version', False)
            if 'module_name' in self.config['tools'][tool]:
                mn = self.config['tools'][tool]['module_name']
                cmd = '%s python load %s' % (self.config['lmod']['path'], mn)
                self.config['tools'][tool].setdefault('module_load', cmd)
                cmd = '%s python unload %s' % (self.config['lmod']['path'], mn)
                self.config['tools'][tool].setdefault('module_unload', cmd)

        uap_tools_path = os.path.join(self._uap_path, 'tools')
        uap_python = os.path.join(self._uap_path, "python_env", "bin", "python")
        for tool_file in os.listdir(uap_tools_path):
            tool_path = os.path.join(uap_tools_path, tool_file)
            if not tool_file.endswith('.py') or not os.path.isfile(tool_path):
                continue
            tool = tool_file[:-3]
            if not tool in self.config['tools'] or \
                    not isinstance(self.config['tools'][tool], dict):
                        self.config['tools'][tool] = dict()
            self.config['tools'][tool].setdefault('path', [uap_python, tool_path])
            self.config['tools'][tool].setdefault('get_version', '--help')
            self.config['tools'][tool].setdefault('exit_code', 0)
            self.config['tools'][tool].setdefault('ignore_version', False)

        if not 'destination_path' in self.config:
            raise UAPError("%s: Missing key: destination_path"
                         % self.config_name)
        if not os.path.exists(self.config['destination_path']):
            raise UAPError("%s: Destination path does not exist: %s"
                         % (self.config_name,
                            self.config['destination_path'])
            )

        # Make self.config['destination_path'] an absolute path if necessary
        if not os.path.isabs(self.config['destination_path']):
            self.config['destination_path'] = os.path.normpath(os.path.join(
                self.get_config_path(), self.config['destination_path']))

        symlink = "%s-out" % self.config['id']
        if not os.path.exists(symlink):
            try:
                os.symlink(self.config['destination_path'], symlink)
            except OSError, e:
                logger.warn(e)
                logger.warn('File %s seems to be a broken symlink and will '
                            'not be overwritten.' % symlink)

        if not 'cluster' in self.config or self.config['cluster'] is None:
            self.config['cluster'] = dict()

        if self.get_cluster_type() is not None:
            self.config['cluster'].setdefault('default_submit_options',
                    self.get_cluster_command('default_options'))
        for i in ['default_submit_options', 'default_pre_job_command',
                'default_post_job_command']:
            self.config['cluster'].setdefault(i, '')
        self.config['cluster'].setdefault('default_job_quota', 0) # no quota

        self.build_steps()

    def get_config(self):
        return self.config

    def build_steps(self):
        self.steps = {}
        if not 'steps' in self.config:
            raise UAPError("%s: Missing key: steps" % self.config_name)
        re_simple_key = re.compile('^[a-zA-Z0-9_]+$')
        re_complex_key = re.compile('^([a-zA-Z0-9_]+)\s+\(([a-zA-Z0-9_]+)\)$')

        # step one: instantiate all steps
        for step_key, step_description in self.config['steps'].items():

            # the step keys in the configuration may be either:
            # - MODULE_NAME 
            # - DIFFERENT_STEP_NAME\s+\(MODULE_NAME\)
            step_name = None
            module_name = None
            if re_simple_key.match(step_key):
                step_name = step_key
                module_name = step_key
            else:
                match = re_complex_key.match(step_key)
                if match:
                    step_name = match.group(1)
                    module_name = match.group(2)

            if step_name == 'temp':
                # A step cannot be named 'temp' because we need the out/temp
                # directory to store temporary files.
                raise UAPError("%s: A step name cannot be 'temp'."
                             % self.config_name)
            step_class = abstract_step.AbstractStep.get_step_class_for_key(module_name)
            step = step_class(self)

            step.set_step_name(step_name)
            step.set_options(step_description)

            self.steps[step_name] = step

        # step two: set dependencies
        for step_name, step in self.steps.items():
            for parent_step in step._options['_depends']:
                if not parent_step in self.steps.keys():
                    raise UAPError("%s: Step %s specifies an undefined "
                                 "dependency: %s."
                                 % (self.config_name,
                                    step_name, parent_step))
                step.add_dependency(self.steps[parent_step])

        # step three: perform topological sort
        # if there's a cycle (yeah, the algorithm is O(n^2), tsk, tsk...)

        unassigned_steps = set(self.steps.keys())
        assigned_steps = set()
        self.topological_step_order = []
        while len(unassigned_steps) > 0:
            # choose all tasks which have all dependencies resolved, either
            # because they have no dependencies or are already assigned
            next_steps = []
            for step_name in unassigned_steps:
                is_ready = True
                for dep in self.steps[step_name].dependencies:
                    dep_name = dep.get_step_name()
                    if not dep_name in assigned_steps:
                        is_ready = False
                        break
                if is_ready:
                    next_steps.append(step_name)
            if len(next_steps) == 0:
                raise UAPError("%s: There is a cycle in the step dependencies."
                             % self.config_name)
            for step_name in misc.natsorted(next_steps):
                self.topological_step_order.append(step_name)
                assigned_steps.add(step_name)
                unassigned_steps.remove(step_name)

        # step four: finalize step
        for step in self.steps.values():
            step.finalize()

    def print_source_runs(self):
        for step_name in self.topological_step_order:
            step = self.steps[step_name]
            if isinstance(step, abstract_step.AbstractSourceStep):
                for run_id in misc.natsorted(step.get_run_ids()):
                    print("%s/%s" % (step, run_id))

    def add_file_dependencies(self, output_path, input_paths):
        if output_path in self.file_dependencies:
            raise UAPError("Different steps/runs/tags want to create "
                         "the same output file: %s." % output_path)
        self.file_dependencies[output_path] = set(input_paths)

        for inpath in input_paths:
            if not inpath in self.file_dependencies_reverse:
                self.file_dependencies_reverse[inpath] = set()
            self.file_dependencies_reverse[inpath].add(output_path)

    def add_task_for_output_file(self, output_path, task_id):
        if output_path in self.task_id_for_output_file:
            raise UAPError("More than one step is trying to create the "
                         "same output file: %s." % output_path)
        self.task_id_for_output_file[output_path] = task_id

        if not task_id in self.output_files_for_task_id:
            self.output_files_for_task_id[task_id] = set()
        self.output_files_for_task_id[task_id].add(output_path)

    def add_task_for_input_file(self, input_path, task_id):
        if not input_path in self.task_ids_for_input_file:
            self.task_ids_for_input_file[input_path] = set()
        self.task_ids_for_input_file[input_path].add(task_id)

        if not task_id in self.input_files_for_task_id:
            self.input_files_for_task_id[task_id] = set()
        self.input_files_for_task_id[task_id].add(input_path)

    def check_command(self, command):
        for argument in command:
            if not isinstance(argument, str):
                raise UAPError(
                    "The command to be launched '%s' contains non-string "
                    "argument '%s'. Therefore the command will fail. Please "
                    "fix this type issue." % (command, argument))
        return

    def setup_lmod(self):
        '''
        If lmod is configured this functions sets the required environmental variables.
        '''
        if 'lmod' in self.config:
            os.environ['MODULEPATH'] = self.config['lmod']['module_path']

    def check_tools(self):
        '''
        checks whether all tools references by the configuration are available
        and records their versions as determined by ``[tool] --version`` etc.
        '''
        if not 'tools' in self.config:
            return
        pool = multiprocessing.Pool()
        def kill_pool(signum, frame):
            pool.close()
            raise UAPError('Keybord interrupt during tool check.')
        original_int_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, kill_pool)
        for tool_id, tool_check_info in \
                pool.imap_unordered(check_tool, self.config['tools'].items()):
            self.tool_versions[tool_id] = tool_check_info
        pool.close()
        signal.signal(signal.SIGINT, original_int_handler)

    def notify(self, message, attachment = None):
        '''
        prints a notification to the screen and optionally delivers the
        message on additional channels (as defined by the configuration)
        '''
        print(message)
        if 'notify' in self.config:
            try:
                notify = self.config['notify']
                match = re.search('^(https?://[^/]+)/([a-z0-9]+)$', notify)
                if match:
                    host = match.group(1)
                    token = match.group(2)
                    args = ['curl', host, '-X', 'POST', '-d', '@-']
                    proc = subprocess.Popen(args, stdin = subprocess.PIPE)
                    data = {'token': token, 'message': message}
                    if attachment:
                        data['attachment_name'] = attachment['name']
                        data['attachment_data'] = base64.b64encode(attachment['data'])
                    proc.stdin.write(json.dumps(data))
                    proc.stdin.close()
                    proc.wait()
                else:
                    logger.warn('Cloud not split patter into http(s)://host/token to notify: %s' %
                            self.config['notify'])
            except:
                # swallow all exception that happen here, failing notifications
                # are no reason to crash the entire thing
                logger.warn('Notification of "%s" failed with:' %
                        (self.config['notify'], sys.exc_info()[0]))
                pass

    def get_cluster_job_ids(self):
        '''
        The argument less method returns a set the cluster job ids of all
        subbmited jobs.
        '''
        ids = set()
        for task in self.all_tasks_topologically_sorted:
            queued_ping_file = task.get_run().get_queued_ping_file()
            failed_qpf = queued_ping_file + '.last'
            if os.path.exists(queued_ping_file):
                with open(queued_ping_file, 'r') as fl:
                    info = yaml.load(fl, Loader=yaml.FullLoader)
                ids.add(info['job_id'])
            elif os.path.exists(failed_qpf):
                with open(failed_qpf, 'r') as fl:
                    info = yaml.load(fl, Loader=yaml.FullLoader)
                ids.add(info['job_id'])
        return ids

    def check_ping_files(self, print_more_warnings = False,
                         print_details = False, fix_problems = False):
        run_problems = list()
        queue_problems = list()
        last_problems = list()
        check_queue = True

        try:
            stat_output = subprocess.check_output(
                [self.get_cluster_command('stat')],
                stderr = subprocess.STDOUT)
        except (KeyError, OSError, subprocess.CalledProcessError):
            # we don't have a stat tool here, if subprocess.CalledProcessError
            # is raised
            check_queue = False

        if print_more_warnings and not check_queue:
            try:
                ce = self.get_cluster_command('stat')
            except KeyError:
                ce = "a cluster engine"
            print("Attention, we cannot check stale queued ping files because "
                  "this host does not have %s." % ce)

        running_jids = set()

        if check_queue:
            for line in stat_output.split("\n"):
                if 'COMPLETING' in line:
                    continue
                try:
                    jid = int(line.strip().split(' ')[0].split('_')[0])
                    running_jids.add(str(jid))
                except ValueError:
                    # this is not a JID
                    pass

        now = datetime.datetime.now()
        for task in self.all_tasks_topologically_sorted:
            exec_ping_file = task.get_run().get_executing_ping_file()
            queued_ping_file = task.get_run().get_queued_ping_file()
            last_queued_ping_file = queued_ping_file + '.last'
            if os.path.exists(exec_ping_file):
                info = yaml.load(open(exec_ping_file, 'r'), Loader=yaml.FullLoader)
                start_time = info['start_time']
                last_activity = datetime.datetime.fromtimestamp(
                    abstract_step.AbstractStep.fsc.getmtime(exec_ping_file))
                last_activity_difference = now - last_activity
                if last_activity_difference.total_seconds() > \
                   abstract_step.AbstractStep.PING_TIMEOUT:
                    run_problems.append((task, exec_ping_file,
                                         last_activity_difference,
                                         last_activity - start_time))

            if os.path.exists(queued_ping_file) and check_queue:
                info = yaml.load(open(queued_ping_file, 'r'), Loader=yaml.FullLoader)
                if not str(info['job_id']) in running_jids:
                    queue_problems.append((task, queued_ping_file,
                                           info['submit_time'], info['job_id']))
            if os.path.exists(last_queued_ping_file):
                info = yaml.load(open(last_queued_ping_file, 'r'), Loader=yaml.FullLoader)
                if not str(info['job_id']) in running_jids:
                    last_problems.append((task, last_queued_ping_file,
                                           info['submit_time'], info['job_id']))

        show_hint = False

        if len(run_problems) > 0:
            show_hint = True
            label = "Warning: There are %d stale run ping files." % len(run_problems)
            print(label)
            if print_details:
                print('-' * len(label))
                run_problems = sorted(run_problems, key=itemgetter(2, 3), reverse=True)
                for problem in run_problems:
                    task = problem[0]
                    path = problem[1]
                    last_activity_difference = problem[2]
                    ran_for = problem[3]
                    print("dead since %13s, ran for %13s: %s" % (
                        misc.duration_to_str(last_activity_difference),
                        misc.duration_to_str(ran_for), task))
                print("")

        if len(queue_problems) > 0:
            show_hint = True
            label = "Warning: There are %d tasks marked as queued, but they do not seem to be queued." % len(queue_problems)
            print(label)
            if print_details:
                print('-' * len(label))
                queue_problems = sorted(queue_problems, key=itemgetter(2), reverse=True)
                for problem in queue_problems:
                    task = problem[0]
                    path = problem[1]
                    start_time = problem[2]
                    job_id = problem[3]
                    print("submitted job %s at %13s: %s" % (job_id, start_time, task))
                print("")

        if len(last_problems) > 0:
            show_hint = True
            label = "Info: %d tasks were queued and ran." % len(last_problems)
            print(label)
            if print_details:
                print('-' * len(label))
                last_problems = sorted(last_problems, key=itemgetter(2), reverse=True)
                for problem in last_problems:
                    task = problem[0]
                    path = problem[1]
                    start_time = problem[2]
                    job_id = problem[3]
                    print("submitted job %s at %13s: %s" % (job_id, start_time, task))
                print("")

        if fix_problems:
            all_problems = run_problems
            all_problems.extend(queue_problems)
            for problem in all_problems:
                path = problem[1]
                print("Now deleting %s..." % path)
                os.unlink(path)

        if show_hint:
            if print_more_warnings and not print_details or not fix_problems:
                print("Hint: Run 'uap %s fix-problems --details' to see the "
                      "details."  % self.args.config.name)
            if not fix_problems:
                print("Hint: Run 'uap %s fix-problems --srsly' to fix these "
                      "problems (that is, delete all problematic ping files)."
                      % self.args.config.name)

    def check_volatile_files(self, details = False, srsly = False):
        collected_files = set()
        for task in self.all_tasks_topologically_sorted:
            collected_files |= task.volatilize_if_possible(srsly)
        if not srsly and len(collected_files) > 0:
            if details:
                for path in sorted(collected_files):
                    print(path)
            total_size = 0
            for path in collected_files:
                total_size += os.path.getsize(path)
            print("Hint: You could save %s of disk space by volatilizing %d "
                  "output files." % (misc.bytes_to_str(total_size),
                                     len(collected_files)))
            print("Call 'uap %s volatilize --srsly' to purge the files."
                  % self.args.config.name)

    def autodetect_cluster_type(self):
        cluster_config = self.get_cluster_config()
        # Let's see if we can successfully run a cluster identity test
        # Test all configured cluster types
        for cluster_type in cluster_config.keys():
            # Do we have an identity test command
            identity = dict()
            for key in ['test', 'answer']:
                try:
                    identity[key] = cluster_config[cluster_type]\
                                        ['identity_%s' % key]
                except KeyError:
                    raise UAPError("%s: Missing 'identity_%s' for %s"
                                 "cluster type."
                                 % (self._cluster_config_path,
                                    key, cluster_type)
                             )
            # Now that we know let's test for that cluster    
            if not isinstance(identity['answer'], list):
                identity['answer'] = [identity['answer']]
            for answer in identity['answer']:
                try:
                    if (subprocess.check_output( identity['test'] )
                        .startswith(answer) ):
                        return cluster_type
                except OSError:
                    pass
        logger.warn('Cluster type could not be detected.')
        return None

    def get_cluster_type(self):
        return self._cluster_type

    def set_cluster_type(self, cluster_type):
        if cluster_type is not None and not cluster_type in self.get_cluster_config():
            raise UAPError ('Cluster type "%s" not configured.' % cluster_type)
        self._cluster_type = cluster_type

    '''
    Shorthand to retrieve a cluster-type-dependent command or filename
    (cc == cluster command).
    '''
    def get_cluster_command(self, key):
        ct = self.get_cluster_type()
        if key not in self.get_cluster_config()[ct].keys():
            raise UAPError('The option "%s" is not available for the cluster "%s".' %
                    (key, ct))
        return self.get_cluster_config()[ct][key]

    '''
    Shorthand to retrieve a cluster-type-dependent command line part (this is a
    list)
    '''
    def get_cluster_command_cli_option(self, key, value):
        result = self.get_cluster_config()[self.get_cluster_type()][key]
        if '%s' in result:
            return [result % value]
        else:
            return [result, value]
