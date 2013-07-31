import sys
from abstract_step import *
import copy
import csv
import glob
import os
import re
import yaml

class FastqSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(FastqSource, self).__init__(pipeline)
        
        self.add_connection('out/reads')

        self.add_option('indices', str, dict, optional = True)
        self.add_option('pattern', str)
        self.add_option('group', str)
        self.add_option('paired_end', bool)
        self.add_option('sample_id_prefix', str, optional = True)
        
    def declare_runs(self):
        regex = re.compile(self.option('group'))
        
        found_files = dict()
        
        # find FASTQ files
        for path in glob.glob(os.path.abspath(self.option('pattern'))):
            match = regex.match(os.path.basename(path))
            if match == None:
                raise StandardError("Couldn't match regex /%s/ to file %s." % (self.option('group'), os.path.basename(path)))
            
            sample_id_parts = []
            if self.option_set_in_config('sample_id_prefix'):
                sample_id_parts.append(self.option('sample_id_prefix'))
                
            sample_id_parts += list(match.groups())
            sample_id = '_'.join(sample_id_parts)
            if not sample_id in found_files:
                found_files[sample_id] = list()
            found_files[sample_id].append(path)

        # declare a run for every sample
        for run_id, paths in found_files.items():
            with self.declare_run(run_id) as run:
                run.add_public_info("paired_end", self.option("paired_end"))
                for path in paths:
                    run.add_output_file("reads", path, [])

        # determine index information...
        # retrieve each run and punch in the information
        if self.option_set_in_config('indices'):
            if type(self.option('indices')) == str:
                # read indices from CSV file
                indices_path = self.option('indices')
                reader = csv.DictReader(open(indices_path))
                for row in reader:
                    sample_id = row['SampleID']
                    run = self.get_run(sample_id)
                    if run != None:
                        index = row['Index']
                        run.add_public_info('index', index)
            else:
                # indices are defined in the configuration
                for sample_id, index in self.option('indices').items():
                    run = self.get_run(sample_id)
                    if run != None:
                        run.add_public_info('index', index)
