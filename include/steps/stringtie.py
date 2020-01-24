from uaperrors import UAPError
import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger = getLogger('uap_logger')


class Stringtie(AbstractStep):

    def __init__(self, pipeline):
        super(Stringtie, self).__init__(pipeline)

        self.set_cores(12)

        self.add_connection('in/alignments')
        self.add_connection('in/reference')
        self.add_connection('out/assembling')
        self.add_connection('out/gene_abund')
        self.add_connection('out/cov_refs')
        self.add_connection('out/log_stderr')
        self.add_connection('out/e2t.ctab')
        self.add_connection('out/e_data.ctab')
        self.add_connection('out/i2t.ctab')
        self.add_connection('out/i_data.ctab')
        self.add_connection('out/t_data.ctab')

        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('stringtie')

        self.add_option('G', str, optional=True, default=None,
                        description="use reference transcript annotation to guide assembly")
        self.add_option('v', bool, optional=True,
                        description='Turns on verbose mode, printing bundle processing details')
        self.add_option('p', int, default=1, optional=True,
                        description='Specify the number of processing threads (CPUs) to use '
                                    'for transcript assembly. The default is 1')
        self.add_option('m', int, optional=True,
                        description='Sets the minimum length allowed for the predicted '
                                    'transcripts. Default: 200')
        self.add_option('l', str, optional=True,
                        description='Sets <label> as the prefix for the name of the output '
                                    'transcripts. Default: STRG')
        self.add_option('f', float, optional=True,
                        description='Sets the minimum isoform abundance of the predicted '
                                    'transcripts as a fraction of the most abundant transcript '
                                    'assembled at a given locus. Lower abundance transcripts are '
                                    'often artifacts of incompletely spliced precursors of '
                                    'processed transcripts. Default: 0.1')

        self.add_option('fr', bool, optional=True,
                        description='assume stranded library fr-secondstrand')

        self.add_option('rf', bool, optional=True,
                        description='assume stranded library fr-firststrand')

        self.add_option('M', float, optional=True,
                        description='Sets the maximum fraction of muliple-location-mapped reads '
                                    'that are allowed to be present at a given locus. Default: '
                                    '0.95.')

        self.add_option('e', bool, optional=True,
                        description=""" Limits the processing of read alignments to only estimate
                                    and output the assembled transcripts matching the reference
                                    transcripts given with the -G option (requires -G, recommended
                                    for -B/-b). With this option, read bundles with no reference
                                    transcripts will be entirely skipped, which may provide a
                                    considerable speed boost when the given set of reference
                                    transcripts is limited to a set of target genes, for
                                    example.""")

        self.add_option('B', bool, optional=True,
                        description=""" This switch enables the output of Ballgown input table
                                    files (*.ctab) containing coverage data for the reference
                                    transcripts given with the -G option. (See the Ballgown
                                    documentation for a description of these files.) With this
                                    option StringTie can be used as a direct replacement of the
                                    tablemaker program included with the Ballgown distribution.
                                    If the option -o is given as a full path to the output
                                    transcript file, StringTie will write the *.ctab files in
                                    the same directory as the output GTF.""")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('p'))

        options = ['v', 'p', 'm', 'l', 'f', 'M', 'e', 'B']

        set_options = [option for option in options if
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                value = str(self.get_option(option))
                if os.path.isfile(value):
                    value = os.path.abspath(value)
                option_list.append('-%s' % option)
                option_list.append(value)

        if self.is_option_set_in_config('fr') and self.get_option('fr'):
            option_list.append('--fr')

        if self.is_option_set_in_config('rf') and self.get_option('rf'):
            option_list.append('--rf')

        # look for reference assembly in in-connections
        ref_assembly = os.path.abspath(self.get_option('G'))
        # check reference annotation
        if ref_assembly is not None and not os.path.isfile(ref_assembly):
            raise UAPError(
                "The path %s provided to option 'G' is not a file."
                % self.get_option('G'))
        for run_id, connection in run_ids_connections_files.items():
            if 'in/reference' in connection.keys():
                if ref_assembly is not None:
                    UAPError('Reference assembly given through option and connection.')
                if len(connection['in/reference']) != 1:
                    UAPError('More then one reference assembly passed from run %s' % run_id)
                ref_assembly = connection['in/reference'][0]
                break


        for run_id, connection in run_ids_connections_files.items():
            # look for run specific reference assembly
            if 'in/alignments' not in connection.keys()\
                    and 'in/reference' in connection.keys():
                continue
            elif 'in/reference' in connection.keys():
                if len(connection['in/reference']) != 1:
                    UAPError('More then one reference assembly passed from run %s' % run_id)
                run_ref_assembly = connection['in/reference'][0]
            else:
                run_ref_assembly = ref_assembly

            alignments = connection['in/alignments']
            # check, if only a single input file is provided
            if len(alignments) != 1:
                raise UAPError("Expected exactly one alignments file %s" % input_paths)
            input_paths = alignments[0]

            run = self.declare_run(run_id)

            assembling = run.add_output_file(
                'assembling',
                '%s-assembling.gtf' % run_id,
                [input_paths])

            log_stderr = run.add_output_file(
                'log_stderr',
                '%s-log_stderr.txt' % run_id,
                [input_paths])

            gene_abund = run.add_output_file(
                'gene_abund',
                '%s-gene_abund.tab' % run_id,
                [input_paths])

            cov_refs = run.add_output_file(
                'cov_refs',
                '%s-cov_refs.gtf' % run_id,
                [input_paths])


            stringtie = [self.get_tool('stringtie'), alignments, '-o', assembling,
                         '-A', gene_abund, '-C', cov_refs]
            if run_ref_assembly is not None:
                stringtie.extend(['-G', run_ref_assembly])
            stringtie.extend(option_list)

            pipe = run.new_exec_group().add_pipeline()
            pipe.add_command(stringtie, stdout_path=assembling, stderr_path=log_stderr)

            if self.is_option_set_in_config('B') and self.get_option('B'):
                mv_exec_group = run.new_exec_group()
                connections = ['e2t.ctab',
                               'e_data.ctab',
                               'i2t.ctab',
                               'i_data.ctab',
                               't_data.ctab']

                for connection in connections:
                    is_produced = ''.join([run.get_output_directory_du_jour_placeholder(),
                                           '/', connection])

                    out_file = run_id + '-' + connection
                    is_wanted = run.add_output_file(connection,
                                                    out_file,
                                                    [input_paths])

                    mv_exec_group.add_command([self.get_tool('mv'),
                                               is_produced, is_wanted])
