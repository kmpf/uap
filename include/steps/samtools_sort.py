from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SamtoolsSort(AbstractStep):
    '''
    automatically recognizes input format
    The step implements samtools sort to sort sam, cram, bam.
    crashes on sam to cram

    Sort alignments by leftmost coordinates, or by read name when -n is used.
    An appropriate @HD-SO sort order header tag will be added or an existing
    one updated if necessary.

    Documentation::

        http://www.htslib.org/doc/samtools.html
    '''

    def __init__(self, pipeline):
        super(SamtoolsSort, self).__init__(pipeline)
        self.set_cores(1)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')

        self.require_tool('dd')
        self.require_tool('samtools')
        # in case of sam output to compress
        self.require_tool('pigz')

        self.add_option('l', str, default  = None,  optional=True,
                        description = 'Set compression level, from 0 (uncompressed) to 9 (best)')
        self.add_option('m', str, default  = None,  optional=True,
                        description = 'Set maximum memory per thread; suffix K/M/G recognized [768M]')
        self.add_option('sort-by-name', bool, default  = None,  optional=False,
                        description = 'sort by read name original option is -n')

        self.add_option('t', str, default  = None,  optional=True,
                        description = 'Sort by value of TAG. Uses position as secondary index (or read name if -n is set)')

        self.add_option('O', str, choices = ['BAM', 'SAM', 'CRAM'], default  = None,  optional=False,
                        description = 'output format bam sam or cram')

        self.add_option('reference', str, default  = None,  optional=True,
                        description = 'reference fasta file need for cram output')


        self.add_option('temp-sort-dir', str, optional = True,
                        description = 'Intermediate sort files are stored intothis directory. original option -T')

        # not samtools specific
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                                    engine and threads ie")

        # [Options for 'dd':]
	self.add_option('fifo', bool, optional=True, default=False,
                description='Enable the FIFO functionality for splitting large input files.')
        self.add_option('dd-blocksize', str, optional=True, default='4096k',
                description='Read data with ``dd`` and set the blocksize.')


    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]

                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise StepError(self, "Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                if self.is_option_set_in_config('temp-sort-dir'):
                    if not os.path.isdir(self.get_option('temp-sort-dir')):
                        #dir not present
                        raise StepError(self, "Directory %s not found" % self.get_option('temp-sort-dir'))
                    if not os.access(self.get_option('temp-sort-dir'), os.W_OK):
                        #not accessible
                        raise StepError(self, "Directory %s not accessible." % self.get_option('temp-sort-dir'))

                    if self.is_option_set_in_config('reference'):
                        if not os.path.isfile(self.get_option('reference')):
                            raise StepError(self, "reference:   %s not accessible" % self.get_option('reference'))

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
                        # 0 use dd to read the data
                        if self.get_option('fifo'):
                            dd_in = [
                                self.get_tool('dd'),
                                'ibs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_paths[0]
                            ]
                            pipe.add_command(dd_in)

                            # 0.1 command: Uncompress file to fifo
                            if is_gzipped:
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--processes', '1',
                                        '--stdout']
                                pipe.add_command(pigz)
                        # 1 command: Sort BAM input
                        samtools_sort = [
                            self.get_tool('samtools'), 'sort',
                            '-O', self.get_option('O')
                        ]

                        if self.is_option_set_in_config('reference'):
                            samtools_sort.extend(['--reference', os.path.abspath(self.get_option('reference'))])


                        if self.get_option('O') == 'CRAM':
                            if  not self.is_option_set_in_config('reference'):
                                raise StepError(self, "option: reference not set in config rqueiered for CRAM output")


                        if self.get_option('sort-by-name'):
                            samtools_sort.append('-n')

                        sortpath =  './'

                        if self.is_option_set_in_config('temp-sort-dir'):
                            sortpath = os.path.join(os.path.abspath(self.get_option('temp-sort-dir')),
                                                'sort' + run_id)

                        samtools_sort.extend(
                            ['-@', str(self.get_option('cores') - 1),
                             '-T', sortpath])


                        suffix = self.get_option('O')
                        postfix = '' if suffix in ['CRAM', 'BAM'] else '.gz'
                        out_path = run.add_output_file(
                            'alignments',
                            '%s.sorted.%s%s'  % (run_id, suffix.lower(), postfix),
                            input_paths)

                        if self.get_option('fifo'):
                            samtools_sort.append('-')
                            pipe.add_command(samtools_sort)
                            samtools_sort = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize')
                            ]
                        else:
                            samtools_sort.extend(input_paths)

                        if suffix in ['CRAM', 'BAM']:
                            # write output directly
                            pipe.add_command(samtools_sort, stdout_path=out_path)
                        else:
                            # output is sam needs to be compressed
                            pipe.add_command(samtools_sort)

                            pigz = [self.get_tool('pigz'),
                                    '--stdout']
                            pipe.add_command(pigz, stdout_path=out_path)

