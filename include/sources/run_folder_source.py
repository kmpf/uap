import sys
from abstract_source import *
import csv
import glob
import os
import pipeline
import yaml

class RunFolderSource(AbstractSource):
    def __init__(self, pipeline, options):
        super(RunFolderSource, self).__init__(pipeline)
        
        raise StandardError("Let's fix this first, shall we?")

        path = options['path']

        self.samples = {}

        if not 'paired_end' in options:
            raise StandardError("missing paired_end key in source")

        if not os.path.exists(path):
            raise ConfigurationException("Source path does not exist: " + path)

        # find all samples
        for sample_path in glob.glob(os.path.join(path, 'Unaligned', 'Project_*', 'Sample_*')):
            sample_name = os.path.basename(sample_path).replace('Sample_', '')
            if sample_name in self.samples:
                raise ConfigurationException("Duplicate sample: " + sample_name)
            self.samples[sample_name] = {}
            self.samples[sample_name]['files'] = []
            for path in sorted(glob.glob(os.path.join(sample_path, '*.fastq.gz'))):
                self.samples[sample_name]['files'].append(path)

            self.samples[sample_name]['info'] = { 'paired_end': options['paired_end'] }

            # read sample sheets
            sample_sheet_path = os.path.join(sample_path, 'SampleSheet.csv')
            reader = csv.DictReader(open(sample_sheet_path))
            for row in reader:
                sample_id = row['SampleID']
                index = row['Index']
                if not 'index' in self.samples[sample_id]['info']:
                    self.samples[sample_id]['info']['index'] = index
                else:
                    if index != self.samples[sample_id]['info']['index']:
                        raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)

            if options['paired_end'] == True:
                # determine R1/R2 info for each input file: read_number
                self.samples[sample_name]['info']['read_number'] = {}
                for path in self.samples[sample_name]['files']:
                    isR1 = '_R1' in path
                    isR2 = '_R2' in path
                    if isR1 and isR2:
                        raise StandardError("Unable to determine read_numer, seems to be both R1 and R2: " + path)
                    if (not isR1) and (not isR2):
                        raise StandardError("Unable to determine read_numer, seems to be neither R1 nor R2: " + path)
                    self.samples[sample_name]['info']['read_number'][os.path.basename(path)] = 'R1' if isR1 else 'R2'
