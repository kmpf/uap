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

        path = options['path']

        self.samples = {}

        # find all samples
        self.all_samples = {}
        if not os.path.exists(path):
            raise ConfigurationException("Source path does not exist: " + path)
        for sample_path in glob.glob(os.path.join(path, 'Unaligned', 'Project_*', 'Sample_*')):
            sample_name = os.path.basename(sample_path).replace('Sample_', '')
            if sample_name in self.all_samples:
                raise ConfigurationException("Duplicate sample: " + sample_name)
            self.samples[sample_name] = {}
            self.samples[sample_name]['files'] = []
            for path in sorted(glob.glob(os.path.join(sample_path, '*.fastq.gz'))):
                self.samples[sample_name]['files'].append(path)

            # read sample sheets
            sample_sheet_path = os.path.join(sample_path, 'SampleSheet.csv')
            reader = csv.DictReader(open(sample_sheet_path))
            for row in reader:
                sample_id = row['SampleID']
                index = row['Index']
                if not 'index' in self.samples[sample_id]:
                    self.samples[sample_id]['index'] = index
                else:
                    if index != self.samples[sample_id]['index']:
                        raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)
