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
                                    files (\*.ctab) containing coverage data for the reference
                                    transcripts given with the -G option. (See the Ballgown
                                    documentation for a description of these files.) With this
                                    option StringTie can be used as a direct replacement of the
                                    tablemaker program included with the Ballgown distribution.
                                    The \*.ctab files will be supplied to child steps through
                                    additional connections ``out/e2t.ctab``, ``out/e_data.ctab``,
                                    ``out/i2t.ctab``, ``out/i_data.ctab`` and ``out/t_data.ctab``.
                                    """)

    def runs(self, cc):
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
        if self.is_option_set_in_config('G'):
            ref_assembly = os.path.abspath(self.get_option('G'))
            if not os.path.isfile(ref_assembly):
                raise UAPError('[Stringtie]: %s is no file.' %
                        self.get_option('G'))
        else:
            ref_assembly = None
        con_ref_assembly = cc.look_for_unique('in/reference', ref_assembly)
        ref_per_run = cc.all_runs_have_connection('in/reference')

        allignment_runs = cc.get_runs_with_connections('in/alignments')
        for run_id in allignment_runs:
            connection = cc[run_id]
            if ref_per_run is True:
                con_ref_assembly = connection['in/reference'][0]

            alignments = connection['in/alignments']
            # check, if only a single input file is provided
            if len(alignments) != 1:
                raise UAPError("Expected exactly one alignments file %s" %
                        alignments)

            run = self.declare_run(run_id)

            assembling = run.add_output_file(
                'assembling',
                '%s-assembling.gtf' % run_id,
                alignments)

            log_stderr = run.add_output_file(
                'log_stderr',
                '%s-log_stderr.txt' % run_id,
                alignments)

            gene_abund = run.add_output_file(
                'gene_abund',
                '%s-gene_abund.tab' % run_id,
                alignments)


            stringtie = [self.get_tool('stringtie'), alignments[0], '-o', assembling,
                         '-A', gene_abund]
            if con_ref_assembly is not None:
                cov_refs = run.add_output_file(
                    'cov_refs',
                    '%s-cov_refs.gtf' % run_id,
                    alignments)
                stringtie.extend(['-C', cov_refs, '-G', ref_assembly])
                if ref_assembly is None:
                    # include dependency
                    alignments.append(con_ref_assembly)
            elif self.is_option_set_in_config('B') \
                    or self.is_option_set_in_config('b') \
                    or self.is_option_set_in_config('e'):
                        UAPError('[stringtie] Options -B, -b and -e can only '
                                'be used if a reference is provided with -G.')
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
                    run.add_out_connection(connection)
                    is_produced = ''.join([run.get_output_directory_du_jour_placeholder(),
                                           '/', connection])

                    out_file = run_id + '-' + connection
                    is_wanted = run.add_output_file(connection,
                                                    out_file,
                                                    alignments)

                    mv_exec_group.add_command([self.get_tool('mv'),
                                               is_produced, is_wanted])
