import copy
import csv
import datetime
import glob
import json
import os
import re
import StringIO
import subprocess
import sys
import yaml

import abstract_step
import abstract_source
import task as task_module

# an enum class, yanked from http://stackoverflow.com/questions/36932/whats-the-best-way-to-implement-an-enum-in-python
class Enum(set):
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

# an exception class for reporting configuration errors
class ConfigurationException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Pipeline(object):

    states = Enum(['WAITING', 'READY', 'FINISHED'])
    run_modes = Enum(['DRY_RUN', 'TEST_RUN', 'FULL'])

    def __init__(self):

        # now determine the Git hash of the repository
        self.git_hash_tag = subprocess.check_output(['git', 'describe', '--all', '--dirty', '--long']).strip()
        if '-dirty' in self.git_hash_tag:
            if not '--even-if-dirty' in sys.argv:
                print("The repository has uncommitted changes, which is why we will exit right now.")
                print("If this is not a production environment, you can skip this test by specifying --even-if-dirty on the command line.")
                exit(1)

        if '--even-if-dirty' in sys.argv:
            sys.argv.remove('--even-if-dirty')
        # the configuration as read from config.yaml
        self.config = {}

        # dictionary of sample names => information
        self.all_samples = {}

        # list of steps, steps are objects with inter-dependencies
        self.steps = []

        # the run mode of the pipeline
        self.run_mode = self.run_modes.FULL
        if '--dry-run' in sys.argv:
            sys.argv.remove('--dry-run')
            self.run_mode = self.run_modes.DRY_RUN
        if '--test-run' in sys.argv:
            sys.argv.remove('--test-run')
            self.run_mode = self.run_modes.TEST_RUN

        self.dry_run_cache = {}

        self.read_config()

        # collect all tasks
        self.task_for_task_id = {}
        self.all_tasks = []
        for step in self.steps:
            for run_id in sorted(step.get_run_ids()):
                task = task_module.Task(self, step, run_id)
                self.all_tasks.append(task)
                self.task_for_task_id[str(task)] = task

        self.tool_versions = {}
        self.check_tools()

        if self.run_mode == self.run_modes.DRY_RUN:
            # fill the dry run cache with input files
            # fill dry_run_cache with source files
            for run_id in self.steps[0].get_run_info().keys():
                for annotation, files in self.steps[0].get_run_info()[run_id]['output_files'].items():
                    for path in files.keys():
                        self.dry_run_cache[path] = datetime.datetime.now()


    # read configuration and make sure it's good
    def read_config(self):
        #print >> sys.stderr, "Reading configuration..."
        self.config = yaml.load(open('config.yaml'))
        if not 'sources' in self.config:
            raise ConfigurationException("Missing key: sources")
        for source in self.config['sources']:
            key = source.keys()[0]
            source_instance = abstract_source.get_source_class_for_key(key)(self, source[key])
            for sample_id, sample_info in source_instance.samples.items():
                if sample_id in self.all_samples:
                    raise ConfigurationException("Sample appears multiple times in sources: " + sample_id)
                self.all_samples[sample_id] = copy.deepcopy(sample_info)
        if not 'destination_path' in self.config:
            raise ConfigurationException("Missing key: destination_path")
        if not os.path.exists(self.config['destination_path']):
            raise ConfigurationException("Destination path does not exist: " + self.config['destination_path'])

        if not os.path.exists("out"):
            os.symlink(self.config['destination_path'], 'out')

        self.build_steps()

    def build_steps(self):
        self.steps = []
        if not 'steps' in self.config:
            raise ConfigurationException("Missing key: steps")

        # stacks are used to handle dependencies between steps
        step_stack = []
        indent_stack = []
        last_indent = None

        source_step = abstract_step.get_step_class_for_key('source')(self)
        self.steps.append(source_step)

        if (self.run_mode == self.run_modes.TEST_RUN):
            # this is a test run, append a head step which boild the
            # problem down to 1000 lines per input file
            head_step = abstract_step.get_step_class_for_key('head')(self)
            head_step.add_dependency(source_step)
            self.steps.append(head_step)
            source_step = head_step
            # now source_step points to head_step and everything else
            # depends on it, HARR HARR!

        steps_definition = self.config['steps']

        # remove comments (from # to EOL)
        while '#' in steps_definition:
            p0 = steps_definition.index('#')
            p1 = steps_definition.index("\n", p0)
            temp = bytearray(steps_definition)
            temp[p0:(p1 + 1)] = ''
            steps_definition = str(temp)

        steps_definition_offset = 0
        while steps_definition_offset < len(steps_definition):
            # find next newline
            newline_index = steps_definition.find("\n", steps_definition_offset)
            # if there's no newline, use the rest of the string
            if newline_index < 0:
                newline_index = len(steps_definition) - 1
            line = steps_definition[steps_definition_offset:(newline_index + 1)]
            steps_definition_offset = newline_index + 1
            if len(line.strip()) == 0:
                continue
            block_mapping = False
            if line.strip()[-1] == '{':
                block_mapping = True
                # line ends with a { --> find matching } and extend line
                level = 1
                while level > 0:
                    c = steps_definition[steps_definition_offset]
                    if c == '{':
                        level += 1
                    elif c == '}':
                        level -= 1
                    line += c
                    steps_definition_offset += 1
            regex = re.compile('(\s*)-\s([^\s]+)(\s\{([^\}]+)\})?')
            match = regex.match(line)
            if match:
                indent = len(match.group(1))
                step_name = match.group(2)
                options_def = match.group(3)
                options = {}
                if options_def:
                    if block_mapping:
                        options_def = options_def.strip()[1:-1]
                    options = yaml.load(options_def)

                # determine class and instatiate it with options
                step_class = abstract_step.get_step_class_for_key(step_name)
                step = step_class(self)
                step.set_options(options)

                # append step to steps list
                self.steps.append(step)

                if last_indent != None:
                    # level_shift > 0 if indent increase, < 0 if decrease, 0 if constant
                    level_shift = indent - last_indent

                    if level_shift > 0:
                        # going down one level, push items on stacks
                        indent_stack.append(last_indent)
                        step_stack.append(self.steps[-2])
                    elif level_shift < 0:
                        # going up, pop items from stack
                        pop_indent = None
                        while len(indent_stack) > 0 and indent_stack[-1] >= indent:
                            pop_indent = indent_stack.pop()
                            pop_step = step_stack.pop()
                        if pop_indent != indent:
                            raise ConfigurationException("invalid indentation: '" + line + "'")

                # if there is an upstream step, add it as a dependency
                # otherwise, link it to the source step which provides the
                # input files
                if len(step_stack) > 0:
                    step.add_dependency(step_stack[-1])
                else:
                    step.add_dependency(source_step)

                last_indent = indent

            else:
                raise ConfigurationException("Invalid steps definition, error in line: '" + line + "'.")

    def print_tasks(self):
        count = {}
        for task in self.all_tasks:
            state = task.get_task_state()
            if not state in count:
                count[state] = 0
            count[state] += 1
            print('[' + task.get_task_state()[0].lower() + '] ' + str(task))
        print('tasks: ' + str(len(self.all_tasks)) + ' total, ' + ', '.join([str(count[_]) + ' ' + _.lower() for _ in sorted(count.keys())]))

    def has_unfinished_tasks(self, task_list):
        unfinished_tasks = [task for task in task_list if task.step.get_run_state(task.run_id) != self.states.FINISHED]
        return len(unfinished_tasks) > 0

    def pick_next_ready_task(self, task_list):
        ready_tasks = [task for task in task_list if task.step.get_run_state(task.run_id) == self.states.READY]
        return ready_tasks[0]

    def check_tools(self):
        if not 'tools' in self.config:
            return
        for tool_id, info in self.config['tools'].items():
            command = [info['path']]
            if 'get_version' in info:
                command.append(info['get_version'])
            exit_code = None
            try:
                proc = subprocess.Popen(command, stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE, close_fds = True)
            except:
                raise ConfigurationException("Tool not found: " + info['path'])
            proc.wait()
            exit_code = proc.returncode
            self.tool_versions[tool_id] = {
                'command': (' '.join(command)).strip(),
                'exit_code': exit_code,
                'response': (proc.stdout.read() + proc.stderr.read()).strip()
            }
            expected_exit_code = 0
            if 'exit_code' in info:
                expected_exit_code = info['exit_code']
            if exit_code != expected_exit_code:
                raise ConfigurationException("Tool check failed for " + tool_id + ": " + ' '.join(command) + ' - exit code is: ' + str(exit_code) + ' (expected ' + str(expected_exit_code) + ')')

    # returns a short description of the configured pipeline
    def __str__(self):
        s = ''
        s += "Number of samples: " + str(len(self.all_samples)) + "\n"
        for sample in sorted(self.all_samples.keys()):
            s += "- " + sample + "\n"
        return s

    def notify(self, message):
        print(message)
        if 'notify' in self.config:
            try:
                notify = self.config['notify']
                match = re.search('^(http://[a-z\.]+:\d+)/([a-z0-9]+)$', notify)
                if match:
                    host = match.group(1)
                    token = match.group(2)
                    args = ['curl', host, '-X', 'POST', '-d', '@-']
                    proc = subprocess.Popen(args, stdin = subprocess.PIPE)
                    proc.stdin.write(json.dumps({'token': token, 'message': message}))
                    proc.stdin.close()
                    proc.wait()
            except:
                # swallow all exception that happen here, failing notifications
                # are no reason to crash the entire thing
                pass
