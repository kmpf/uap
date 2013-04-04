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

        if not 'paired_end' in options:
            raise StandardError("missing paired_end key in source")

        for path in glob.glob(pattern):
            match = regex.match(os.path.basename(path))
            sample_id_parts = []
            if 'sample_id_prefix' in options:
                sample_id_parts.append(options['sample_id_prefix'])
            sample_id_parts += list(match.groups())
            sample_id = '_'.join(sample_id_parts)
            if not sample_id in self.samples:
                self.samples[sample_id] = { 'files': [], 'info': { 'paired_end': options['paired_end'] } }
            self.samples[sample_id]['files'].append(path)

        if options['indices'].__class__ == str:
            # read indices from CSV file or dictionary
            indices_path = options['indices']
            reader = csv.DictReader(open(indices_path))
            for row in reader:
                sample_id = row['SampleID']
                if sample_id in self.samples:
                    index = row['Index']
                    if not 'index' in self.samples[sample_id]['info']:
                        self.samples[sample_id]['info']['index'] = index
                    else:
                        if index != self.samples[sample_id]['info']['index']:
                            raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)
        else:
            # indices are defined in the configuration
            for sample_id, index in options['indices'].items():
                if sample_id in self.samples:
                    if not 'index' in self.samples[sample_id]['info']:
                        self.samples[sample_id]['info']['index'] = index
                    else:
                        if index != self.samples[sample_id]['info']['index']:
                            raise StandardError("Inconsistent index defined in sample sheets for sample " + sample_id)

        if options['paired_end'] == True:
            # determine R1/R2 info for each input file: read_number
            for sample_name in self.samples.keys():
                self.samples[sample_name]['info']['read_number'] = {}
                for path in self.samples[sample_name]['files']:
                    isR1 = '_R1' in path
                    isR2 = '_R2' in path
                    if isR1 and isR2:
                        raise StandardError("Unable to determine read_numer, seems to be both R1 and R2: " + path)
                    if (not isR1) and (not isR2):
                        raise StandardError("Unable to determine read_numer, seems to be neither R1 nor R2: " + path)
                    self.samples[sample_name]['info']['read_number'][os.path.basename(path)] = 'R1' if isR1 else 'R2'
