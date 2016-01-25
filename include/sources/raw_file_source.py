import sys
from abstract_step import *
import os
import glob

class RawFileSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(RawFileSource, self).__init__(pipeline)

        self.add_connection('out/raw')

        self.add_option(
            'pattern', str, optional = True,
            description = "A file name pattern, for example "
            "``/home/test/fastq/Sample_*.fastq.gz``.")

        self.add_option(
            'group', str, optional = True, 
            description = "A regular expression which is applied to found "
            "files, and which is used to determine the sample name from the "
            "file name. For example, `(Sample_\d+)_R[12].fastq.gz``, when "
            "applied to a file called ``Sample_1_R1.fastq.gz``, would result "
            "in a sample name of ``Sample_1``. You can specify multiple "
            "capture groups in the regular expression.")

        self.add_option('sample_id_prefix', str, optional = True )

        self.add_option('sample_to_files_map', dict, str, 
                        description = "A listing of sample names and their "
                        "associated files. This must be provided as a YAML "
                        "dictionary.", optional = True)

    def runs(self, run_ids_connections_files):
        # found_files holds the runIDs and their related files
        found_files = dict()

        if self.is_option_set_in_config('group') and \
           self.is_option_set_in_config('pattern'):
            regex = re.compile(self.get_option('group'))
            
            # find files matching the 'group' pattern in all files matching 'pattern'
            for path in glob.glob(os.path.abspath(self.get_option('pattern'))):
                match = regex.match(os.path.basename(path))
                if match == None:
                    raise StandardError("Couldn't match regex /%s/ to file %s."
                        % (self.get_option('group'), os.path.basename(path)))
            
                sample_id_parts = []
                if self.is_option_set_in_config('sample_id_prefix'):
                    sample_id_parts.append(self.get_option('sample_id_prefix'))
                
                sample_id_parts += list(match.groups())
                sample_id = '_'.join(sample_id_parts)
                if not sample_id in found_files:
                    found_files[sample_id] = list()
                found_files[sample_id].append(path)
                
        elif self.is_option_set_in_config('sample_to_files_map'):
            for run_id, paths in self.get_option('sample_to_files_map'):
                for path in paths:
                    if not os.path.isfile(path):
                        raise StandardError("[raw_file_source]: %s is no file. "
                                            "Please provide correct path." 
                                            % path)
                if not run_id in found_files:
                    found_files[run_id] = list()
                found_files[run_id] = paths

        else:
            raise StandardError("[raw_file_source]: Either 'group' AND 'pattern'"
                                " OR 'sample_to_files_map' options have to be "
                                "set. ")

        # declare a run for every sample
        for run_id, paths in found_files.items():
            with self.declare_run(run_id) as run:
                for path in paths:
                    run.add_output_file("raw", path, [])
