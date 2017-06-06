from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class AdapterRemoval(AbstractStep):
    '''
    AdapterRemoval (ver. 2.1.7)

    This program searches for and removes remnant adapter sequences from
    your read data.  The program can analyze both single end and paired end
    data.  For detailed explanation of the parameters, please refer to the
    man page.  For comments, suggestions  and feedback please contact Stinus
    Lindgreen (stinus@binf.ku.dk) and Mikkel Schubert (MikkelSch@gmail.com).

    If you use the program, please cite the paper:
    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid
    adapter trimming, identification, and read merging.
    BMC Research Notes, 12;9(1):88.

    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2

    "Pipeline specific "input and output expected to be gzipped"
    '''

    def __init__(self, pipeline):
        super(AdapterRemoval, self).__init__(pipeline)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        self.add_connection('out/collapsed')
        self.add_connection('out/collapsed.truncated')
        self.add_connection('out/discarded')
        self.add_connection('out/truncated')
        self.add_connection('out/pair1.truncated')
        self.add_connection('out/pair2.truncated')
        self.add_connection('out/singleton.truncated')
        self.add_connection('out/settings')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.require_tool('adapterremoval')
        self.require_tool('pwd')
        self.require_tool('mv')

        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                    engine and threads ie")

        self.add_option('treatAs', str, optional=False, default=None,
                        choices=['paired', 'single'],
                        description="Pipeline specific: Sets how input is \
                        fastq files are treated: paired R1 and R2 in one \
                        command, single r1 and r2 get clipped independently \
                        if none r2 just empty")

        self.add_option('qualitybase', str, optional=True, default=None,
                        choices=['33', '64', 'solexa'],
                        description="Quality base used to encode Phred scores \
                        in input; either 33, 64, or solexa [current: 33]")

        self.add_option('adapter1', str, optional=False,
                        description="Adapter sequence expected to be found in \
                        mate 1 reads [current: \
                        AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG]")

        self.add_option('adapter2', str, optional=False,
                        description="Adapter sequence expected to be found in \
                        mate 2 reads [current: \
                        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT]")

        self.add_option('mm', int, default=None, optional=True,
                        description="Max error-rate when aligning reads \
                        and/or adapters. If > 1, the max error-rate is set \
                        to 1 / MISMATCH_RATE; if < 0, the defaults are used, \
                        otherwise the user-supplied value is used directly. \
                        [defaults: 1/3 for trimming; 1/10 when identifing \
                        adapters]")

        self.add_option('maxns', int, default=None, optional=True,
                        description="Reads containing more ambiguous bases \
                        (N) than this number after trimming are discarded \
                        [current: 1000]")

        self.add_option('shift', int, default=None, optional=True,
                        description="Consider alignments where up to N \
                        nucleotides are missing from the 5' termini \
                        [current: 2]")

        self.add_option('trimns', bool, default=None, optional=True,
                        description="If set, trim ambiguous bases (N) \
                        at 5'/3' termini [current: off]")

        self.add_option('trimqualities', bool, default=None, optional=True,
                        description="If set, trim bases at 5'/3' termini \
                        with quality scores <= to --minquality value \
                        [current: off]")

        self.add_option('minquality', int, default=None, optional=True,
                        description="Inclusive minimum; see --trimqualities \
                        for details [current: 2]")

        self.add_option('minlength', int, default=None, optional=True,
                        description="Reads shorter than this length are \
                        discarded following trimming [current: 15]")

        self.add_option('maxlength', int, default=None, optional=True,
                        description="Reads longer than this length are \
                        discarded following trimming [current: 4294967295]")

        self.add_option('collapse', bool, default=None, optional=False,
                        description="When set, paired ended read alignments \
                        of --minalignmentlength or more bases are combined \
                        into a single consensus sequence, representing the \
                        complete insert, and written to either \
                        basename.collapsed or basename.collapsed.truncated \
                        (if trimmed due to low-quality bases following \
                        collapse); for single-ended reads, putative complete \
                        inserts are identified as having at least \
                        --minalignmentlength bases overlap with the adapter \
                        sequence, and are written to the the same files \
                        [current:off]")

        self.add_option('minalignmentlength', bool, default=None,
                        optional=True, description="If --collapse is set, \
                        paired reads must overlap at least this number of \
                        bases to be collapsed, and single-ended reads must \
                        overlap at least this number of bases with the \
                        adapter to be considered complete template molecules \
                        [current: 11]")

        self.add_option('minadapteroverlap', bool, default=None, optional=True,
                        description="In single-end mode, reads are only \
                        trimmed if the overlap between read and the adapter \
                        is at least X bases long, not counting ambiguous \
                        nucleotides (N); this is independant of the \
                        --minalignmentlength when using --collapse, allowing \
                        a conservative selection of putative complete inserts \
                        while ensuring that all possible adapter \
                        contamination is trimmed [current: 0].")

        self.add_option('identify-adapters', bool, default=None, optional=True)
        self.add_option('seed', int, default=22595, optional=True)
        self.add_option('threads', int, default=1, optional=True)

    def runs(self, run_ids_connections_files):

        self.__treat_as_paired = True if \
            self.get_option('treatAs') == 'paired' \
            else False

        self.set_cores(self.get_option('cores'))
        print(run_ids_connections_files)
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # set everything empty
                out_connections = ['collapsed', 'collapsed.truncated',
                                   'discarded', 'settings']

                if self.__treat_as_paired:
                    out_connections.extend([
                        'pair1.truncated',
                        'pair2.truncated',
                        'singleton.truncated'])
                else:
                    out_connections.append('truncated')

                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                if self.__treat_as_paired:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]

                ar_exec_group = run.new_exec_group()
                ar = [self.get_tool('adapterremoval')]
                ar.extend(['--file1', r1, ])

                if self.__treat_as_paired:
                    ar.extend(['--file2', r2, ])

                ar.extend(['--adapter1', self.get_option('adapter1')])

                if self.__treat_as_paired:
                    ar.extend(['--adapter2', self.get_option('adapter2')])

                basename = run.get_output_directory_du_jour_placeholder() + '/' + run_id
                ar.extend(['--basename', basename])

                ar.extend(['--gzip'])

                if self.is_option_set_in_config('threads'):
                    ar.extend(['--threads', str(self.get_option('threads'))])

                if self.is_option_set_in_config('qualitybase'):
                    ar.extend(['--qualitybase',
                               str(self.get_option('qualitybase'))])

                if self.get_option('collapse') is True:
                    ar.extend(['--collapse'])

                if self.is_option_set_in_config('mm'):
                    ar.extend(['--mm',  str(self.get_option('mm'))])

                if self.is_option_set_in_config('maxns'):
                    ar.extend(['--maxns',  str(self.get_option('maxns'))])

                if self.is_option_set_in_config('shift'):
                    ar.extend(['--shift',  str(self.get_option('shift'))])

                if self.is_option_set_in_config('trimns'):
                    ar.extend(['--trimns'])

                if self.is_option_set_in_config('trimqualities'):
                    ar.extend(['--trimqualities'])

                if self.is_option_set_in_config('minquality'):
                    ar.extend(['--minquality',
                               str(self.get_option('minquality'))])

                if self.is_option_set_in_config('minlength'):
                    ar.extend(['--minlength',
                               str(self.get_option('minlength'))])

                if self.is_option_set_in_config('maxlength'):
                    ar.extend(['--maxlength',
                               str(self.get_option('maxlength'))])

                if self.is_option_set_in_config('minalignmentlength'):
                    ar.extend(['--minalignmentlength',
                               str(self.get_option('minalignmentlength'))])

                if self.is_option_set_in_config('minadapteroverlap'):
                    ar.extend(['--minadapteroverlap',
                               str(self.get_option('minadapteroverlap'))])

                if self.is_option_set_in_config('seed'):
                    ar.extend(['--seed', str(self.get_option('seed'))])

                output_fileset = [r1]
                if self.__treat_as_paired:
                    output_fileset.append(r2)

                stderr_file = "%s-adapterremoval-log_stderr.txt" % (run_id)
                log_stderr = run.add_output_file("log_stderr",
                                                 stderr_file, output_fileset)
                stdout_file = "%s-adapterremoval-log_stdout.txt" % (run_id)
                log_stdout = run.add_output_file("log_stdout",
                                                 stdout_file, output_fileset)

                mv_exec_group = run.new_exec_group()

                for connection in out_connections:
                    if self.get_option('collapse') is False:
                        if connection == 'collapsed' or \
                           connection == 'collapsed.truncated':
                            continue

                    if connection == 'settings':
                        settings_file = run_id + '.' + connection
                        run.add_output_file(connection,
                                            settings_file, output_fileset)
                    else:
                        out_file = run_id + '.' + connection + '.fastq.gz'
                        iswanted = run.add_output_file(connection,
                                                       out_file,
                                                       output_fileset)
                        isproduced = run.get_output_directory_du_jour_placeholder() + '/' + run_id + '.' + connection + '.gz'

                        mv_exec_group.add_command([self.get_tool('mv'),
                                                   isproduced, iswanted])

                ar_exec_group.add_command(ar, stdout_path=log_stdout,
                                          stderr_path=log_stderr)
