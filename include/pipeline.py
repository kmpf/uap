import csv
import glob
import os
import re
import StringIO
import sys
import yaml

import abstract_step

# an exception class for reporting configuration errors
class ConfigurationException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Pipeline(object):

    def __init__(self):

        # the configuration as read from config.yaml
        self.config = {}

        # dictionary of sample names => information
        self.all_samples = {}

        # list of steps
        self.steps = []

        self.read_config()
        self.gather_information()

    # read configuration and make sure it's good
    def read_config(self):
        print >> sys.stderr, "Reading configuration..."
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
        self.build_steps()

    # for every source path, look for samples in [path]/Unaligned/Project_*/Sample_*
    def gather_information(self):
        print >> sys.stderr, "Gathering information..."

        # find all samples
        self.all_samples = {}
        for path in self.config['sourcePaths']:
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

        for line in StringIO.StringIO(self.config['steps']):
            #sys.stdout.write(line)
            regex = re.compile('(\s*)-\s([^\s]+)(\s\{([^\}]+)\})?')
            match = regex.match(line)
            if match:
                indent = len(match.group(1))
                step_name = match.group(2)
                options = yaml.load(match.group(3)) if match.group(3) else {}
                print(indent, step_name, options)

                # determine class and instatiate it with options
                step_class = abstract_step.get_step_class_for_key(step_name)
                step = step_class()
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
                if len(step_stack) > 0:
                    step.add_dependency(step_stack[-1])

                last_indent = indent

            else:
                raise ConfigurationException("Invalid steps definition, error in line: '" + line + "'.")
        for step in self.steps:
            print(abstract_step.AbstractStep.__str__(step))
        #for key in self.config['steps']:
            #self.steps.append(abstract_step.get_step_class_for_key(key))

    # returns a short description of the configured pipeline
    def __str__(self):
        s = ''
        s += "Number of samples: " + str(len(self.all_samples)) + "\n"
        for sample in sorted(self.all_samples.keys()):
            s += "- " + sample + "\n"
        return s
