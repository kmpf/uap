import sys
from abstract_step import *
import csv
import glob
import string
import os
import pipeline
import yaml

class RunFolderSource(AbstractSourceStep):
    '''
    This source looks for fastq.gz files in ``[path]/Unaligned/Project_*/Sample_*`` 
    and pulls additional information from CSV sample sheets it finds. It also makes 
    sure that index information for all samples is coherent and unambiguous.
    '''
    
    def __init__(self, pipeline):
        super(RunFolderSource, self).__init__(pipeline)
        
        self.add_connection('out/reads')
        
        self.add_option('path', str)
        self.add_option('paired_end', bool)
        
    def declare_runs(self):
        path = self.get_option('path')

        if not self.is_option_set_in_config('paired_end'):
            raise StandardError("missing paired_end key in source")

        if not os.path.exists(path):
            raise StandardError("Source path does not exist: " + path)

        found_samples = dict()

        # find all samples
        for sample_path in glob.glob(os.path.join(path, 'Unaligned', 'Project_*', 'Sample_*')):
            sample_name = os.path.basename(sample_path).replace('Sample_', '')
            if sample_name in found_samples:
                raise StandardError("Duplicate sample: " + sample_name)
            
            if not sample_name in found_samples:
                found_samples[sample_name] = list()

            for path in sorted(glob.glob(os.path.join(sample_path, '*.fastq.gz'))):
                found_samples[sample_name].append(path)

        # create a run for each sample in found_samples and store data therein
        for run_id, paths in found_samples.items():
            with self.declare_run(run_id) as run:
                run.add_public_info('paired_end', self.get_option('paired_end'))
                for path in paths:
                    run.add_output_file('reads', path, [])
                    
                if self.get_option('paired_end') == True:
                    # determine R1/R2 info for each input file: read_number
                    r1_files = list()
                    r2_files = list()
                    for path in paths:
                        isR1 = '_R1' in path
                        isR2 = '_R2' in path
                        if isR1 and isR2:
                            raise StandardError("Unable to determine read_number, " +
                                                "seems to be both R1 and R2: " + path)
                        if (not isR1) and (not isR2):
                            raise StandardError("Unable to determine read_number, " +
                                                "seems to be neither R1 nor R2: " + path)
                        if isR1:
                            r1_files.append(os.path.basename(path))
                        if isR2:
                            r2_files.append(os.path.basename(path))

                    run.add_public_info('R1', r1_files)
                    run.add_public_info('R2', r2_files)



                # read sample sheets
                sample_path =  (os.path.dirname(paths[0]))
                sample_sheet_path = os.path.join(sample_path, 'SampleSheet.csv')
                reader = csv.DictReader(open(sample_sheet_path))
                # get and set indices
                for row in reader:
                    sample_id = row['SampleID']
                    if not sample_id in found_samples.keys():
                        raise StandardError("Found sample %s in %s, but it shouldn't be here." 
                                            % sample_id, sample_sheet_path)

                    index = row['Index'].split('-')

                    if len(index) == 2:
                        if not run.has_public_info('index-R1') and not run.has_public_info('index-R2'):
                            run.add_public_info('index-R1', index[0])
                            run.add_public_info('index-R2', index[1])
                        else: 
                            stored_index = (run.get_public_info('index-R1') + '-' 
                                            + run.get_public_info('index-R2'))
                            if index.join('-') != stored_index:
                                raise StandardError("Inconsistent index defined in %s for sample %s"
                                                % (sample_sheet_path, sample_id))

                    elif len(index) == 1:
                        if not run.has_public_info('index-R1'):
                            run.add_public_info('index-R1', index[0])
                        elif index[0] != run.get_public_info('index-R1'):
                            raise StandardError("Inconsistent index defined in %s for sample %s"
                                                % (sample_sheet_path, sample_id))

                    else:
                        raise StandardError("Unknown index definition %s found in %s" % 
                                            (index.join('-'), sample_sheet_path))
