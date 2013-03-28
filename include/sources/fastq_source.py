import sys
from abstract_source import *
import csv
import glob
import os
import pipeline
import re
import yaml

class FastqSource(AbstractSource):
    def __init__(self, pipeline, options):
        super(FastqSource, self).__init__(pipeline)

        pattern = options['pattern']
        regex = re.compile(options['group'])

        self.samples = {}

        for path in glob.glob(pattern):
            match = regex.match(os.path.basename(path))
            sample_id_parts = []
            if 'sample_id_prefix' in options:
                sample_id_parts.append(options['sample_id_prefix'])
            sample_id_parts += list(match.groups())
            sample_id = '_'.join(sample_id_parts)
            if not sample_id in self.samples:
                self.samples[sample_id] = { 'files': [] }
            self.samples[sample_id]['files'].append(path)

        indices_path = options['indices']
        reader = csv.DictReader(open(indices_path))
        for row in reader:
            sample_id = row['SampleID']
            index = row['Index']
            if not 'index' in self.samples[sample_id]:
                self.samples[sample_id]['index'] = index
            else:
                if index != self.samples[sample_id]['index']:
                    raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)
