from uaperrors import StepError
import sys
from abstract_step import *
import copy
import csv
import glob
import os
import re
import yaml

import misc


class FastqSource(AbstractSourceStep):

    '''
    The FastqSource class acts as a source for FASTQ files. This source creates a
    run for every sample.

    Specify a file name pattern in *pattern* and define how sample names should
    be determined from file names by specifyign a regular expression in *group*.

    Sample index barcodes may specified by providing a filename to a CSV file
    containing the columns *Sample_ID* and *Index* or directly by defining a
    dictionary which maps indices to sample names.
    '''

    def __init__(self, pipeline):
        super(FastqSource, self).__init__(pipeline)

        self.add_connection('out/first_read')
        self.add_connection('out/second_read', optional=True)

        self.add_option('pattern', str, optional=True,
                        description="A file name pattern, for example "
                        "``/home/test/fastq/Sample_*.fastq.gz``.")

        self.add_option(
            'group',
            str,
            optional=True,
            description="A regular expression which is applied to found "
            "files, and which is used to determine the sample name from "
            r"the file name. For example, ``(Sample_\d+)_R[12].fastq.gz``, "
            "when applied to a file called ``Sample_1_R1.fastq.gz``, would "
            "result in a sample name of ``Sample_1``. You can specify "
            "multiple capture groups in the regular expression.")

        self.add_option('paired_end', bool, description="Specify whether "
                        "the samples are paired end or not.", optional=True)

        self.add_option(
            'indices',
            str,
            dict,
            optional=True,
            description="path to a CSV file or a dictionary of sample_id: "
            "barcode entries.")

        self.add_option(
            'sample_id_prefix',
            str,
            optional=True,
            description="This optional prefix is prepended to every sample "
            "name.")

        self.add_option('sample_to_files_map', dict, str, optional=True,
                        description="A listing of sample names and their "
                        "associated files. This must be provided as a YAML "
                        "dictionary.")

        self.add_option('first_read', str, optional=False,
                        description="Part of the file name that marks all "
                        "files containing sequencing data of the first read. "
                        "Example: 'R1.fastq' or '_1.fastq'")

        self.add_option('second_read', str, optional=True, default=None,
                        description="Part of the file name that marks all "
                        "files containing sequencing data of the second read. "
                        "Example: 'R2.fastq' or '_2.fastq'")

    def runs(self, run_ids_connections_files):
        # found_files holds the runIDs and their related files
        found_files = dict()
        read_types = {'first_read': self.get_option('first_read')}

        paired_end = self.is_option_set_in_config('second_read')
        if self.is_option_set_in_config('paired_end'):
            if paired_end and not self.get_option('paired_end'):
                raise StepError(self,
                    'Second read passed but paired_end set to False.')
            elif not paired_end and self.get_option('paired_end'):
                raise StepError(self,
                    'No second read passed but paired_end set to True.')
        if paired_end:
            read_types['second_read'] = self.get_option('second_read')

        if self.is_option_set_in_config(
                'group') and self.is_option_set_in_config('pattern'):
            regex = re.compile(self.get_option('group'))

            # find FASTQ files
            for path in glob.glob(os.path.abspath(self.get_option('pattern'))):
                match = regex.match(os.path.basename(path))
                if match is None:
                    raise Exception(
                        "Couldn't match regex /%s/ to file %s." %
                        (self.get_option('group'), os.path.basename(path)))

                # Construct sample_id
                sample_id_parts = []
                if self.is_option_set_in_config('sample_id_prefix'):
                    sample_id_parts.append(self.get_option('sample_id_prefix'))

                sample_id_parts += list(match.groups())
                sample_id = '_'.join(sample_id_parts)

                if sample_id not in found_files:
                    found_files[sample_id] = dict()
                # either finds a single match or throws an error
                which_read = misc.assign_string(os.path.basename(path),
                                                read_types.values())
                # so this is save because which_read can only be a value from
                # read_types
                if which_read not in found_files[sample_id]:
                    found_files[sample_id][which_read] = list()
                found_files[sample_id][which_read].append(path)

        elif self.is_option_set_in_config('sample_to_files_map'):
            for sample_id, paths in self.get_option(
                    'sample_to_files_map').items():
                for path in paths:
                    if not os.path.isfile(path):
                        raise Exception("[fastq_source]: %s is no file. "
                                        "Please provide correct path."
                                        % path)

                    if sample_id not in found_files:
                        found_files[sample_id] = dict()
                        # either finds a single match or throws an error
                    which_read = misc.assign_string(os.path.basename(path),
                                                    read_types.values())
                    # so this is save because which_read can only be a value
                    # from read_types
                    if which_read not in found_files[sample_id]:
                        found_files[sample_id][which_read] = list()
                    found_files[sample_id][which_read].append(path)

        else:
            raise Exception("[raw_file_source]: Either 'group' AND 'pattern'"
                            " OR 'sample_to_files_map' options have to be "
                            "set. ")

        # declare a run for every sample
        for run_id in found_files.keys():
            with self.declare_run(run_id) as run:
                run.add_public_info("paired_end", paired_end)
                for read in ['first_read', 'second_read']:
                    if read in read_types.keys():
                        for path in found_files[run_id][read_types[read]]:
                            run.add_output_file(read, path, [])

                # save public information
                run.add_public_info(
                    "first_read", self.get_option("first_read"))
                if paired_end:
                    run.add_public_info(
                        "second_read", self.get_option("second_read"))

        # determine index information...
        # retrieve each run and punch in the information
        if self.is_option_set_in_config('indices'):
            if isinstance(self.get_option('indices'), str):
                # read indices from CSV file
                indices_path = self.get_option('indices')
                reader = csv.DictReader(open(indices_path))
                for row in reader:
                    sample_id = row['SampleID']
                    run = self.get_run(sample_id)
                    if run is not None:
                        index = row['Index'].split('-')
                        if not (run.has_public_info('index-R1')
                                or run.has_public_info('index-R2')):
                            if len(index) == 2:
                                run.add_public_info('index-R1', index[0])
                                run.add_public_info('index-R2', index[1])
                            elif len(index) == 1:
                                run.add_public_info('index-R1', index[0])
                            else:
                                raise Exception(
                                    "Index %s is not a valid index in %s" %
                                    index.join('-'), indices_path)
            else:
                # indices are defined in the configuration
                for sample_id, index in self.get_option('indices').items():
                    run = self.get_run(sample_id)
                    if run is not None:
                        idx = index.split('-')
                        if not (run.has_public_info('index-R1')
                                or run.has_public_info('index-R2')):
                            if len(idx) == 2:
                                run.add_public_info('index-R1', idx[0])
                                run.add_public_info('index-R2', idx[1])
                            elif len(idx) == 1:
                                run.add_public_info('index-R1', idx[0])
                            else:
                                raise Exception(
                                    "Index %s is not a valid index in %s" %
                                    idx.join('-'), index)
