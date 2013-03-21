import copy
import csv
import datetime
import glob
import os
import re
import StringIO
import subprocess
import sys
import yaml

import abstract_step
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

    def __init__(self, run_mode = run_modes.FULL):

        # the configuration as read from config.yaml
        self.config = {}

        # dictionary of sample names => information
        self.all_samples = {}

        # list of steps, steps are objects with inter-dependencies
        self.steps = []

        # the run mode of the pipeline
        self.run_mode = run_mode

        self.read_config()

        # collect all tasks
        self.all_tasks = []
        for step in self.steps:
            for run_id in step.get_run_ids():
                task = task_module.Task(self, step, run_id)
                self.all_tasks.append(task)

        self.check_tools()

    # read configuration and make sure it's good
    def read_config(self):
        #print >> sys.stderr, "Reading configuration..."
        self.config = yaml.load(open('config.yaml'))
        if not 'sourcePaths' in self.config:
            raise ConfigurationException("Missing key: sourcePaths")
        for path in self.config['sourcePaths']:
            if not os.path.exists(path):
                raise ConfigurationException("Source path does not exist: " + path)
        if not 'destinationPath' in self.config:
            raise ConfigurationException("Missing key: destinationPath")
        if not os.path.exists(self.config['destinationPath']):
            raise ConfigurationException("Destination path does not exist: " + self.config['destinationPath'])

        self.gather_information()
        self.build_steps()

    # for every source path, look for samples in [path]/Unaligned/Project_*/Sample_*
    def gather_information(self):
        #print >> sys.stderr, "Gathering information..."

        # find all samples
        self.all_samples = {}
        for path in self.config['sourcePaths']:
            if not os.path.exists(path):
                raise ConfigurationException("Source path does not exist: " + path)
            for samplePath in glob.glob(os.path.join(path, 'Unaligned', 'Project_*', 'Sample_*')):
                sample_name = os.path.basename(samplePath).replace('Sample_', '')
                if sample_name in self.all_samples:
                    raise ConfigurationException("Duplicate sample: " + sample_name)
                self.all_samples[sample_name] = {}
                self.all_samples[sample_name]['path'] = samplePath

        # read sample sheets
        for sample_name, sample in self.all_samples.items():
            sample_sheet_path = os.path.join(sample['path'], 'SampleSheet.csv')
            reader = csv.DictReader(open(sample_sheet_path))
            self.all_samples[sample_name]['lanes'] = {}
            for row in reader:
                self.all_samples[sample_name]['lanes'][row['Lane']] = row

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

                if step_name == 'source':
                    raise ConfigurationException("You cannot use 'source' as a step, it's included automatically.")

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
        print("task states: [w]aiting, [r]eady, [f]inished")
        for task in self.all_tasks:
            print(task)

    def dry_run(self):
        temp_task_list = copy.copy(self.all_tasks)

        # the dry_run_cache contains all files which are created during the process
        dry_run_cache = {}
        # fill dry_run_cache with source files
        for run_id, files in self.steps[0].get_run_info().items():
            for path in files.keys():
                dry_run_cache[path] = datetime.datetime.now()

        while len(temp_task_list) > 0:
            # find a task which is ready
            ready_tasks = [_ for _ in temp_task_list if _['step'].get_run_state(_['run_id'], dry_run_cache) == 'r']
            if len(ready_tasks) == 0:
                raise StandardError("Unable to find a task which is ready.")
            picked_task = ready_tasks[0]
            # now dry-run the taks
            print("now dry running: " + str(picked_task))
            temp_task_list = [_ for _ in temp_task_list if _ != picked_task]
            picked_task['step'].dry_run(picked_task['run_id'], dry_run_cache)
        print("Dry run finished successfully.")

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
            command = [info['path'], info['get_version']]
            exit_code = None
            with open(os.devnull, 'w') as devnull:
                exit_code = subprocess.call(command, stdout = devnull, stderr = devnull)
            if exit_code != 0:
                raise ConfigurationException("Tool check failed for " + tool_id + ": " + ' '.join(command))

    # returns a short description of the configured pipeline
    def __str__(self):
        s = ''
        s += "Number of samples: " + str(len(self.all_samples)) + "\n"
        for sample in sorted(self.all_samples.keys()):
            s += "- " + sample + "\n"
        return s
