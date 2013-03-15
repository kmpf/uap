import csv
import glob
import os
import sys
import yaml

class Pipeline:

    # an exception class for reporting configuration errors
    class ConfigurationException(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)
        
    def __init__(self):
        
        # the configuration as read from config.yaml
        self.config = {}
        
        # dictionary of sample names => information
        self.all_samples = {}
        
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

    def __str__(self):
        s = ''
        s += 'Pipeline information\n'
        s += '--------------------\n'
        s += "Number of samples: " + str(len(self.all_samples)) + ".\n"
        for sample in sorted(self.all_samples.keys()):
            s += sample + "\n"
        return s
    