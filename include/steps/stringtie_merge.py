import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class StringTieMerge(AbstractStep):

    '''StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. merge is a mode of the StringTie tool that is used to assemble transcripts from multiple input files (assemblies). It generates a unified non-redundant set of isoforms.

    NOTE: This step implements the merging part of stringtie. If you want
    stringtie to assemble transcripts from multiple BAM files please use step stringtie!

    https://ccb.jhu.edu/software/stringtie/

    '''
    def __init__(self, pipeline):
        super(StringTieMerge, self).__init__(pipeline)

        self.set_cores(6)

        # The transcripts that should be merged
        self.add_connection('in/features',
                            constraints ={'min_files_per_run': 1, 'max_files_per_run': 1})
        # A .gtf file used as guide for assembly
        self.add_connection('in/features',
                            constraints = {'total_files': 1})

        self.add_connection('out/features')   # contains the assempled transcripts (GTF), -o
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.require_tool('stringtie')

        ## options for stringtie program
        # -G <FILE.gtf/gff>
        self.add_option('G', str, optional = True,
                        description = 'reference annotation to use for guiding the assembly process '
                        '(GTF/GFF3)')
        # -m <min_len>
        self.add_option('m', int, optional = True,
                        description = 'minimum input transcript length to include in the merge (default: 50)')
        # -c <min_cov>
        self.add_option('c', float, optional = True,
                        description = 'minimum input transcript coverage to include in the merge (default: 0)')
        # -F <min_fpkm>
        self.add_option('F', float, optional = True,
                        description = 'minimum input transcript FPKM to include in the merge (default: 1.0)')
        # -T <min_tpm>
        self.add_option('T', float, optional = True,
                        description = 'minimum input transcript TPM to include in the merge (default: 1.0)')
        # -f <min_iso>
        self.add_option('f', float, optional = True,
                        description = 'minimum isoform fraction (default: 0.01)')
        # -g <gap_len>
        self.add_option('g', int, optional = True,
                        description = 'gap between transcripts to merge together (default: 250)')
        # -i
        self.add_option('i', bool, optional = True,
                        description = 'keep merged transcripts with retained introns; by default')
        # -l <label>
        self.add_option('l', str, optional=True,
                        description='name prefix for output transcripts (default: MSTRG)')

        self.add_option('run_id', str, optional=True, default="mergeGTF",
                        description="A name for the run. Since this step merges multiple assemblies (gtf files) "
                        "into a single one, the run_id cannot be the sample name anymore.")

    def runs(self, run_ids_connections_files):

        # Compile the list of options
        options=['m', 'c', 'F', 'T', 'f', 'g', 'i', 'l']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append( '-%s' % option )
                option_list.append( str(self.get_option(option)) )

        features_path = str

        if self.is_option_set_in_config('G'):
            features_path = self.get_option('G')
        else:
            for run_id in run_ids_connections_files.keys():
                try:
                    features_path = run_ids_connections_files[run_id]['in/features'][0]
                except KeyError:
                    continue

        if features_path:
            option_list.append('-G')
            option_list.append(features_path)

        merge_id = self.get_option('run_id')

        input_list = list()

        for run_id in run_ids_connections_files.keys():

            input_paths = run_ids_connections_files[run_id]['in/features']
            input_list.append(input_paths[0])

        with self.declare_run(merge_id) as run:
            outfile = run.add_output_file('features',
                                          '%s-merged-transcripts.gtf' % merge_id,
                                          input_paths)
            stdout = run.add_output_file('log_stdout',
                                         '%s-stringtie_merge.stdout' % merge_id,
                                         input_paths)
            stderr = run.add_output_file('log_stderr',
                                         '%s-stringtie_merge.stderr' % merge_id,
                                         input_paths)

            with run.new_exec_group() as exec_group:

                stringtie_merge = [self.get_tool('stringtie'),
                                   '--merge']

                stringtie_merge.extend(option_list)

                stringtie_merge.extend(['-o', outfile])

                stringtie_merge.extend(input_list)

                exec_group.add_command(stringtie_merge,
                                       stdout_path = stdout,
                                       stderr_path = stderr)
