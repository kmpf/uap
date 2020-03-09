from uaperrors import StepError
import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger = getLogger('uap_logger')


class StringtieMerge(AbstractStep):

    '''
    # stringtie --merge <gtf.list> > outputpat/outputname

    StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. merge is a mode of the StringTie tool that is used to assemble transcripts from multiple input files (assemblies). It generates a unified non-redundant set of isoforms.

    NOTE: This step implements the merging part of stringtie. If you want
    stringtie to assemble transcripts from multiple BAM files please use step stringtie!

    https://ccb.jhu.edu/software/stringtie/
    '''

    def __init__(self, pipeline):
        super(StringtieMerge, self).__init__(pipeline)

        self.set_cores(2)

        # all .gft assemblies from all samples that have been produced with
        # stringtie
        self.add_connection('in/features', format=['gtf', 'gff3'],
                            description='Feature annotations to be merged.')
        self.add_connection(
            'in/reference',
            format=[
                'gtf',
                'gff3'],
            optional=True,
            description='Reference assembly. Can also be passed with option G '
            'or left out for denovo assembling.')
        # merged assembly 'merged.gft'
        self.add_connection('out/features', format='gtf')  # merged.gtf
        self.add_connection('out/assemblies')  # input assemblies txt file
        self.add_connection('out/log_stderr')

        self.require_tool('stringtie')
        self.require_tool('printf')
        self.require_tool('mkdir')
        self.require_tool('mv')

        self.add_option(
            'G',
            str,
            optional=True,
            description='reference annotation to include in the merging (GTF/GFF3)')
        self.add_option(
            'm',
            int,
            optional=True,
            description='minimum input transcript length to include in the merge (default: 50)')
        self.add_option(
            'c',
            int,
            optional=True,
            description='minimum input transcript coverage to include in the merge (default: 0)')
        self.add_option(
            'F',
            float,
            optional=True,
            description='minimum input transcript FPKM to include in the merge (default: 1.0)')
        self.add_option(
            'T',
            float,
            optional=True,
            description='minimum input transcript TPM to include in the merge (default: 1.0)')
        self.add_option('f', float, optional=True,
                        description='minimum isoform fraction (default: 0.01)')
        self.add_option(
            'g',
            int,
            optional=True,
            description='gap between transcripts to merge together (default: 250)')
        self.add_option(
            'i',
            bool,
            optional=True,
            description='keep merged transcripts with retained introns; by default')
        self.add_option(
            'l',
            str,
            optional=True,
            description='name prefix for output transcripts (default: MSTRG)')

        self.add_option('p', int, optional=True,
                        default=2, description='Number of cores')

        self.add_option(
            'output_prefix',
            str,
            optional=True,
            default="merge",
            description='Prefix used in the utput directory.')

    def runs(self, cc):

        # reset cores to number of threads
        self.set_cores(self.get_option('p'))

        # compile list of options
        options = ['m', 'c', 'F', 'T', 'f', 'g', 'i', 'l']

        set_options = [option for option in options if
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                value = str(self.get_option(option))
                option_list.append('-%s' % option)
                option_list.append(value)

        ref_assembly = self.get_option('G')
        if ref_assembly is not None:
            ref_assembly = os.path.abspath(ref_assembly)
        con_ref_assembly = cc.look_for_unique('in/reference', ref_assembly)
        if cc.all_runs_have_connection('in/reference'):
            raise StepError(
                self, 'For stringtieMerge only one reference assmbly can be used.')

        input_files = []
        if con_ref_assembly is not None:
            option_list.extend(['-G', con_ref_assembly])
            if ref_assembly is None:
                # include dependency
                input_files.append(con_ref_assembly)

        # get all paths to the stringtie assemblies from each sample

        stringtie_sample_gtf = []
        assembling_runs = cc.get_runs_with_connections('in/features')
        for run_id in assembling_runs:
            stringtie_sample_gtf.append(cc[run_id]['in/features'][0])

        run_id = self.get_option('output_prefix')
        run = self.declare_run(run_id)

        # create the filename of the assemblies.txt file
        assemblies = [self.get_tool('printf'), '\n'.join(stringtie_sample_gtf)]
        # print assemblies

        input_files.extend(stringtie_sample_gtf)
        assemblies_file = run.add_output_file(
            'assemblies', '%s-stringtieMerge-assemblies.txt' %
            run_id, input_files)

        # 1. create assemblies file
        with run.new_exec_group() as exec_group:
            exec_group.add_command(assemblies, stdout_path=assemblies_file)
            with exec_group.add_pipeline() as stringtie_pipe:
                res = run.add_output_file('features',
                                          '%s-stringtieMerge-merged.gtf' %
                                          run_id, input_files)

                log_err_file = run.add_output_file(
                    'log_stderr', '%s-stringtieMerge-log_stderr.txt' %
                    run_id, input_files)

                stringtieMerge = [self.get_tool('stringtie'), '--merge']
                stringtieMerge.extend(option_list)
                stringtieMerge.append(assemblies_file)
                stringtie_pipe.add_command(stringtieMerge,
                                           stderr_path=log_err_file,
                                           stdout_path=res)
