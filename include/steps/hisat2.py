import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Hisat2(AbstractStep):
    '''
    HISAT2 is a fast and sensitive alignment program for mapping
    next-generation sequencing reads (both DNA and RNA) to a population of
    human genomes (as well as to a single reference genome).

    https://ccb.jhu.edu/software/hisat2/index.shtml
    '''

    def __init__(self, pipeline):
        super(Hisat2, self).__init__(pipeline)
        self.set_cores(12)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/log_stderr')
        self.add_connection('out/metrics')
        self.add_connection('out/unaligned')

        self.require_tool('pigz')
        self.require_tool('hisat2')

        self.add_option('index', str, optional=False,
                        description="Path to hisat2 index (not containing \
                        file suffixes).")

        self.add_option('cores', int, default=12)

        flags = ["q", "qseq", "f", "c", "ignore-quals", "nofw", "dta",
                 "norc", "no-mixed", "no-discordant", "quiet",
                 "qc-filter",
                 "non-deterministic", "remove-chrname", "add-chrname",
                 "no-temp-splicesite", "no-spliced-alignment", "tmo",
                 "no-head", "no-sq", "omit-sec-seq", 'un-gz']

        for flag in flags:
            self.add_option(flag, bool, default=None, optional=True)

        # processing everything as strings is easier
        strflags = ["n-ceil", "ma", "mp", "sp", "np", "rdg",
                    "rdf", "score-min", "k", "rg",
                    "pen-cansplice", "pen-noncansplice", "pen-canintronlen",
                    "pen-noncanintronlen", "min-intronlen", "max-intronlen",
                    "known-splicesite-infile", "novel-splicesite-outfile",
                    "novel-splicesite-infile"]

        for flag in strflags:
            self.add_option(flag, str, default=None, optional=True)

        self.add_option('fr', bool, default=None, optional=False)
        self.add_option('rf', bool, default=None, optional=False)
        self.add_option('ff', bool, default=None, optional=False)

        self.add_option('rg-id', bool, default=None, optional=True,
                        description="puts sample name in rg")

        # note to self just allow R F U and if paired in extend accordingly
        # otherwise mixed single paired will not work
        # Truseq is RF (R)
        self.add_option('rna-strandness', str, choices=["R", "F", "U"],
                        default=None, optional=False,
                        desription="paried and are extended F-> FR, R->RF")

    def runs(self, run_ids_connections_files):
        flags = ["q", "qseq", "f", "c", "ignore-quals", "nofw", "dta",
                 "norc", "no-mixed",  "no-discordant", "quiet", "qc-filter",
                 "non-deterministic", "remove-chrname", "add-chrname",
                 "no-temp-splicesite", "no-spliced-alignment", "tmo",
                 "no-head", "no-sq", "omit-sec-seq"]

        strflags = ["n-ceil", "ma", "mp", "sp", "np", "rdg",
                    "rdf", "score-min", "k", "rg",
                    "pen-cansplice", "pen-noncansplice", "pen-canintronlen",
                    "pen-noncanintronlen", "min-intronlen", "max-intronlen",
                    "known-splicesite-infile", "novel-splicesite-outfile",
                    "novel-splicesite-infile"]

        self.set_cores(self.get_option('cores'))

        res = [self.get_option('fr'),
               self.get_option('rf'),
               self.get_option('ff')]

        if sum(res) > 1:
            message = "too many stranded flags fr, rf, ff: %s"
            raise StandardError(message % (res))

        # Check if option values are valid
        if not os.path.exists(self.get_option('index') + '.1.ht2'):
            logger.error("Could not find index file: %s.*" %
                         self.get_option('index'))
            sys.exit(1)

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read'][0]
                sr_input = run_ids_connections_files[run_id]['in/second_read'][0]

                is_paired_end = True
                input_paths = [fr_input]

                if sr_input is None:
                    is_paired_end = False
                else:
                    input_paths.append(sr_input)

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as hisat2_pipe:
                        # Assemble hisat2 command
                        hisat2 = [self.get_tool('hisat2')]

                        for flag in flags:
                            if self.is_option_set_in_config(flag):
                                if self.get_option(flag) is True:
                                    if flag in ['q', 'f', 'r', 'c']:
                                        hisat2.extend(['-' + flag])
                                    else:
                                        hisat2.extend(['--' + flag])

                        for flag in strflags:
                            if self.is_option_set_in_config(flag):
                                hisat2.extend(['--' + flag,
                                              self.get_option(flag)])

                        hisat2.extend(['-p', str(self.get_option('cores') - 2),
                                       '-x', self.get_option('index')])

                        if is_paired_end:
                            if self.get_option('rna-strandness') == 'F':
                                hisat2.extend(['--rna-strandness', 'FR'])
                            elif self.get_option('rna-strandness') == 'R':
                                hisat2.extend(['--rna-strandness', 'RF'])

                            hisat2.extend([
                                '-1', fr_input,
                                '-2', sr_input])
                        else:
                            hisat2.extend(['-U', fr_input])
                            if not self.get_option('rna-strandness') == 'U':
                                hisat2.extend([
                                    '--rna-strandness',
                                    self.get_option('rna-strandness')])

                        log_stderr = run.add_output_file(
                            'log_stderr',
                            '%s-hisat2-log_stderr.txt' % run_id,
                            input_paths)

                        metrics = run.add_output_file(
                            'metrics',
                            '%s-hisat2-metrics.txt' % run_id,
                            input_paths)
                        hisat2.extend(['--met-file', metrics])

                        if self.is_option_set_in_config('un-gz'):
                            if self.get_option('un-gz') is True:
                                unaligned = run.add_output_file(
                                    'unaligned',
                                    '%s-hisat2-unaligned.fastq.gz' % run_id,
                                    input_paths)
                                hisat2.extend(['--un-gz', unaligned])

                        hisat2_pipe.add_command(hisat2, stderr_path=log_stderr)
                        res = run.add_output_file(
                            'alignments',
                            '%s-hisat2-results.sam.gz' % run_id,
                            input_paths)

                        # Compress hisat2 output
                        pigz = [self.get_tool('pigz'),
                                '--stdout']
                        hisat2_pipe.add_command(pigz, stdout_path=res)
