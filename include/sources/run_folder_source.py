import sys
import csv
import glob
import logging
import string
import os
import yaml

from abstract_step import *
import pipeline

logger = logging.getLogger("uap_logger")

class RunFolderSource(AbstractSourceStep):
    '''
    This source looks for fastq.gz files in 
    ``[path]/Unaligned/Project_*/Sample_*`` and pulls additional information
    from CSV sample sheets it finds. It also makes sure that index information
    for all samples is coherent and unambiguous.
    '''
    
    def __init__(self, pipeline):
        super(RunFolderSource, self).__init__(pipeline)
        
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        
        self.add_option('path', str)
        self.add_option('project', str, default='*')


        self.add_option('paired_end', bool)
        self.add_option('first_read', str, default = "_R1",
            description = "Part of the file name that marks all "
                "files containing sequencing data of the first read. "
                "Example: '_R1.fastq' or '_1.fastq'")

        self.add_option('second_read', str, default = "_R2",
            description = "Part of the file name that marks all "
                "files containing sequencing data of the second read. "
                "Example: 'R2.fastq' or '_2.fastq'")
        
    def runs(self, run_ids_connections_files):

        found_samples = dict()
        read_types = dict()
        if self.get_option('first_read'):
            read_types['first_read'] = self.get_option('first_read')

        if self.get_option('second_read'):
            read_types['second_read'] = self.get_option('second_read')

        # let's look if the path to our data exists 
        path = self.get_option('path')
        path = os.path.abspath(path)
        if not os.path.exists(path):
            raise StandardError("Source path does not exist: " + path)

        # find all samples
        
        project = 'Project_' + self.get_option('project') 
        for sample_path in glob.glob(os.path.join(path, 'Unaligned', project, 'Sample_*')):
            
            sample_name = os.path.basename(sample_path).replace('Sample_', '')
            if sample_name in found_samples:
                raise StandardError("Duplicate sample: " + sample_name)
            
            if not sample_name in found_samples:
                found_samples[sample_name] = dict()

            for path in sorted(glob.glob(os.path.join(sample_path, '*.fastq.gz'))):
                which_read = misc.assign_string(os.path.basename(path),
                                                read_types.values())

                if not which_read in found_samples[sample_name]:
                    found_samples[sample_name][which_read] = list()
                found_samples[sample_name][which_read].append(path)

        # create a run for each sample in found_samples and store data therein
        for run_id in found_samples.keys():
            with self.declare_run(run_id) as run:
                run.add_public_info('paired_end', self.get_option('paired_end'))
                
                sample_path = None
                for read in ['first_read', 'second_read']:
                    if read in read_types.keys():
                        for path in found_samples[run_id][read_types[read]]:
                            run.add_output_file(read, path, [])
                            sample_path = os.path.dirname(path)
                            
                    # always set the out connection even for zero files
                    else:
                        run.add_empty_output_connection(read)

                # read sample sheets
                sample_sheet_path = os.path.join(sample_path, 'SampleSheet.csv')
                csv_file = open(sample_sheet_path)
                reader = csv.DictReader(csv_file)

                # get and set indices
                for row in reader:
                    sample_id = row['SampleID']
                    if not sample_id in found_samples.keys():
                        raise StandardError("Found sample %s in %s, but it "
                                            "shouldn't be here." 
                                            % sample_id, sample_sheet_path)

                    index = row['Index'].split('-')

                    if len(index) == 2:
                        if not run.has_public_info('index-R1') and not run.has_public_info('index-R2'):
                            run.add_public_info('index-R1', index[0])
                            run.add_public_info('index-R2', index[1])
                        else: 
                            stored_index = (run.get_public_info('index-R1') + '-' 
                                            + run.get_public_info('index-R2'))

                            if ('-').join(index)                            != stored_index:
                                raise StandardError("Inconsistent index defined "
                                    "in %s for sample %s" % (sample_sheet_path,
                                                             sample_id))

                    elif len(index) == 1:
                        if not run.has_public_info('index-R1'):
                            run.add_public_info('index-R1', index[0])
                        elif index[0] != run.get_public_info('index-R1'):
                            raise StandardError("Inconsistent index defined in "
                                "%s for sample %s" % (sample_sheet_path, 
                                                      sample_id))

                    else:
                        raise StandardError("Unknown index definition %s found "
                            "in %s" % (('-').join(index), sample_sheet_path))

                csv_file.close()
