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


class Stringtie(AbstractStep):

    '''
    StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step
    to assemble and quantitate full-length transcripts representing multiple splice variants for
    each gene locus. Its input can include not only the alignments of raw reads used by other
    transcript assemblers, but also alignments longer sequences that have been assembled from those
    reads.In order to identify differentially expressed genes between experiments, StringTie's
    output can be processed by specialized software like Ballgown, Cuffdiff or other programs
    (DESeq2, edgeR, etc.).

    NOTE: This step implements that part of stringtie that assembles new transcripts. If you want
    stringtie to assemble transcripts from multiple input files please use step stringtie_merge!

    https://ccb.jhu.edu/software/stringtie/

    '''

    def __init__(self, pipeline):
        super(Stringtie, self).__init__(pipeline)

        self.set_cores(6)

        self.add_connection('in/alignments', format='bam')
        self.add_connection(
            'in/features',
            format=[
                'gtf',
                'gff3'],
            optional=True,
            description='Reference assembly. Can also be passed with option G '
            'or left out for denovo assembling.')
        self.add_connection(
            'out/features',
            format='gtf',
            description='Contains the assempled transcripts (-o).')
        self.add_connection('out/abundances', format='tab',
                            description='Feature abundancies (-A).')
        self.add_connection(
            'out/coverage',
            optional=True,
            format='gtf',
            description='Coverage of the reference assmbly (-B, requires -G)')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        self.add_connection(
            'out/e_data',
            optional=True,
            format='ctab',
            description='Ballgown output (requires -G and -B).')
        self.add_connection(
            'out/i2t',
            optional=True,
            format='ctab',
            description='Ballgown output (requires -G and -B).')
        self.add_connection(
            'out/e2t',
            optional=True,
            format='ctab',
            description='Ballgown output (requires -G and -B).')
        self.add_connection(
            'out/i_data',
            optional=True,
            format='ctab',
            description='Ballgown output (requires -G and -B).')
        self.add_connection(
            'out/t_data',
            optional=True,
            format='ctab',
            description='Ballgown output (requires -G and -B).')

        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('mkfifo')
        self.require_tool('dd')
        self.require_tool('stringtie')

        # -G <FILE.gtf/gff>
        self.add_option(
            'G',
            str,
            optional=True,
            default=None,
            description='reference annotation to use for guiding the assembly process')
        # -v
        self.add_option(
            'v',
            bool,
            optional=True,
            description='Turns on verbose mode, printing bundle processing details')
        # -p <int>
        self.add_option(
            'p',
            int,
            default=6,
            optional=True,
            description='Specify the number of processing threads (CPUs) to use '
            'for transcript assembly.')
        # -m <int>
        self.add_option(
            'm',
            int,
            optional=True,
            description='Sets the minimum length allowed for the predicted '
            'transcripts. Default: 200')
        # -a <INT>
        self.add_option(
            'a',
            int,
            optional=True,
            description='minimum anchor length for junctions (default: 10)')
        # -j <FLOAT>
        self.add_option('j', float, optional=True,
                        description='minimum junction coverage (default: 1)')
        # -t
        self.add_option(
            't',
            bool,
            optional=True,
            description='disable trimming of predicted transcripts based on coverage '
            '(default: coverage trimming is enabled)')
        # -c <FLOAT>
        self.add_option(
            'c',
            float,
            optional=True,
            description='minimum reads per bp coverage to consider for transcript '
            'assembly (default: 2.5)')
        # -g <INT>
        self.add_option(
            'g',
            int,
            optional=True,
            description='gap between read mappings triggering a new bundle (default: 50)')
        # -l <label>
        self.add_option(
            'l',
            str,
            optional=True,
            description='Sets <label> as the prefix for the name of the output '
            'transcripts. Default: STRG')
        # -f <0.1-1.0>
        self.add_option(
            'f',
            float,
            optional=True,
            description='Sets the minimum isoform abundance of the predicted '
            'transcripts as a fraction of the most abundant transcript '
            'assembled at a given locus. Lower abundance transcripts are '
            'often artifacts of incompletely spliced precursors of '
            'processed transcripts. Default: 0.1')
        # --fr
        self.add_option('fr', bool, optional=True,
                        description='assume stranded library fr-secondstrand')
        # --rf
        self.add_option('rf', bool, optional=True,
                        description='assume stranded library fr-firststrand')
        # -M <0.0-1.0>
        self.add_option(
            'M',
            float,
            optional=True,
            description='Sets the maximum fraction of muliple-location-mapped reads '
            'that are allowed to be present at a given locus. Default: '
            '0.95.')
        # -e
        self.add_option(
            'e',
            bool,
            optional=True,
            description='Limits the processing of read alignments to only estimate '
            'and output the assembled transcripts matching the reference '
            'transcripts given with the -G option (requires -G, recommended '
            'for -B/-b). With this option, read bundles with no reference '
            'transcripts will be entirely skipped, which may provide a '
            'considerable speed boost when the given set of reference '
            'transcripts is limited to a set of target genes, for '
            'example.')
        # -B
        self.add_option(
            'B',
            bool,
            optional=True,
            description='This switch enables the output of Ballgown input table '
            r'files (\*.ctab) containing coverage data for the reference '
            'transcripts given with the -G option. (See the Ballgown '
            'documentation for a description of these files.) With this '
            'option StringTie can be used as a direct replacement of the '
            'tablemaker program included with the Ballgown distribution. '
            r'The \*.ctab files will be supplied to child steps through '
            'additional connections ``out/e2t``, ``out/e_data``, '
            '``out/i2t``, ``out/i_data`` and ``out/t_data``.')
        # -b
        self.add_option('b', bool, optional=True,
                        description='enable output of Ballgown table files '
                        'but these files will be created under '
                        'the directory path given as <dir_path>')
        # -x <seqid_list>
        self.add_option(
            'x',
            str,
            optional=True,
            description='Ignore all read alignments (and thus do not attempt to '
            'perform transcript assembly) on the specified reference sequences. '
            'Parameter <seqid_list> can be a single reference sequence name (e.g. '
            '-x chrM) or a comma-delimited list of sequence names (e.g. -x '
            '"chrM,chrX,chrY"). This can speed up StringTie especially in the case '
            'of excluding the mitochondrial genome, whose genes may have very high '
            'coverage in some cases, even though they may be of no interest for a '
            'particular RNA-Seq analysis. The reference sequence names are case '
            'sensitive, they must match identically the names of chromosomes/contigs '
            'of the target genome against which the RNA-Seq reads were aligned in '
            'the first place.')

        # [Options for 'dd':]
        self.add_option(
            'fifo',
            bool,
            optional=True,
            default=False,
            description='Enable the FIFO functionality for splitting large input files.')
        self.add_option('dd-blocksize', str, optional=True, default="2M",
                        description='Provide the blocksize for dd tool.')

    def runs(self, cc):
        self.set_cores(self.get_option('p'))

        options = ['v', 'p', 'm', 'l', 'f', 'M', 'e', 'B']
        options = [
            'l',
            'f',
            'm',
            'a',
            'j',
            't',
            'c',
            'v',
            'g',
            'M',
            'p',
            'B',
            'e',
            'x']

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

        if self.get_option('fr') is not None:
            option_list.append('--fr')

        if self.get_option('rf') is not None:
            option_list.append('--rf')

        # look for reference assembly in in-connections
        option_ref_assembly = self.get_option('G')
        if option_ref_assembly is not None:
            option_ref_assembly = os.path.abspath(option_ref_assembly)
            if not os.path.isfile(option_ref_assembly):
                raise StepError(self, '%s is no file.' %
                                self.get_option('G'))
        ref_assembly = cc.look_for_unique('in/features', option_ref_assembly)
        ref_per_run = cc.all_runs_have_connection('in/features')

        alignment_runs = cc.get_runs_with_connections('in/alignments')
        for run_id in alignment_runs:
            connection = cc[run_id]
            if ref_per_run:
                ref_assembly = connection['in/features'][0]

            alignments = connection['in/alignments']
            # check, if only a single input file is provided
            if len(alignments) != 1:
                raise StepError(
                    self, "Expected exactly one alignments file %s" %
                    alignments)

            run = self.declare_run(run_id)

            assembling = run.add_output_file(
                'features',
                '%s-assembling.gtf' % run_id,
                alignments)

            log_stderr = run.add_output_file(
                'log_stderr',
                '%s-log_stderr.txt' % run_id,
                alignments)

            log_stdout = run.add_output_file(
                'log_stdout',
                '%s-log_stdout.txt' % run_id,
                alignments)

            abundances = run.add_output_file(
                'abundances',
                '%s-abundances.tab' % run_id,
                alignments)

            exec_group = run.new_exec_group()
            if self.get_option('fifo'):
                # 1. create FIFO for BAM file
                fifo_path_bam = run.add_temporary_file('bam_path_fifo',
                                                       designation='input')
                mkfifo_bam = [self.get_tool('mkfifo'), fifo_path_bam]
                exec_group.add_command(mkfifo_bam)

                # 2. read BAM and output to FIFO
                dd_bam = [self.get_tool('dd'),
                          'bs=%s' % self.get_option('dd-blocksize'),
                          'if=%s' % alignments[0],
                          'of=%s' % fifo_path_bam]
                exec_group.add_command(dd_bam)

                # 3. initialize the stringtie command on FIFO
                stringtie = [self.get_tool('stringtie'), fifo_path_bam]
            else:
                # 1. initialize the stringtie command on input BAM file
                stringtie = [self.get_tool('stringtie'), alignments[0]]

            stringtie.extend(['-o', assembling, '-A', abundances])
            if ref_assembly is not None:
                cov_refs = run.add_output_file(
                    'coverage',
                    '%s-cov_refs.gtf' % run_id,
                    alignments)
                stringtie.extend(['-C', cov_refs, '-G', option_ref_assembly])
                if option_ref_assembly is None:
                    # include dependency
                    alignments.append(ref_assembly)
            elif self.is_option_set_in_config('B') \
                    or self.is_option_set_in_config('b') \
                    or self.is_option_set_in_config('e'):
                raise StepError(self, '[stringtie] Options -B, -b and -e can '
                                'only be used if a reference is '
                                'provided with -G.')
            stringtie.extend(option_list)

            exec_group.add_command(stringtie, stdout_path=log_stdout,
                                   stderr_path=log_stderr)

            if self.get_option('B') is not None:
                # rename the ballgown output
                connections = ['e2t',
                               'e_data',
                               'i2t',
                               'i_data',
                               't_data']

                for connection in connections:
                    is_produced = connection + '.ctab'

                    out_file = run_id + '-' + connection + '.ctab'
                    is_wanted = run.add_output_file(connection,
                                                    out_file,
                                                    alignments)

                    run.new_exec_group().add_command([self.get_tool('mv'),
                                                      is_produced, is_wanted])
