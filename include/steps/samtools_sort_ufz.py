from uaperrors import UAPError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SamtoolsSortUFZ(AbstractStep):
    '''
    The step samtools_sort builds on 'samtools sort' to sort SAM files and
    output SAM/BAM/CRAM files. 

    Sort alignments by leftmost coordinates, or by read name when -n is used.
    
    An appropriate @HD-SO sort order header tag will be added or an existing
    one updated if necessary.

    Documentation::

        http://www.htslib.org/doc/samtools.html
    '''

    def __init__(self, pipeline):
        super(SamtoolsSort, self).__init__(pipeline)
        
        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('pigz')

        # [Options for 'samtools':]
        self.add_option('compression-level', int, optional = True, default = 0,
                        description = 'Set compression level, from 0 (uncompressed) '
                        'to 9 (best)')
        self.add_option('max-mem-per-thread', str, optional = True, default = '768M',
                        description = 'Set maximum memory per thread; '
                        'suffix K/M/G recognized [768M]')
        self.add_option('sort-by-name', bool, default = False, optional = True,
                        description = '')
        self.add_option('output-format', str, optional = False, default = 'sam',
                        description = 'Write output as FORMAT (sam/bam/cram)')
        self.add_option('genome-faidx', str, optional = False, 
                        description = 'samtools index of the corresponding genome')
        self.add_option('temp-sort-dir', str, optional = False,
                        description = 'Intermediate sort files are stored into'
                        'this directory.')

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "4096k")

    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise UAPError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                if self.is_option_set_in_config('temp-sort-dir'):
                    if not os.path.isdir(self.get_option('temp-sort-dir')):
                        #dir not present
                        raise UAPError("Directory %s not found" % self.get_option('temp-sort-dir'))
                    if not os.access(self.get_option('temp-sort-dir'), os.W_OK):
                        #not accessible
                        raise UAPError("Directory %s not accessible." % self.get_option('temp-sort-dir'))
                

                    with run.new_exec_group() as exec_group:

                        with exec_group.add_pipeline() as pipe:
                            # 1. command: Read file in chunks
                            dd_in = [
                                self.get_tool('dd'),
                                'ibs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_paths[0]
                            ]
                            pipe.add_command(dd_in)

                            # 1.1 command: Uncompress file to fifo
                            if is_gzipped:
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--processes', '1',
                                        '--stdout']
                                pipe.add_command(pigz)

                            # 2. sort the sam file
                     
                            samtools_sort = [
                                self.get_tool('samtools'), 'sort',
                                '-O', self.get_option('output-format')
                            ]

                            if self.get_option('sort-by-name'):
                                samtools_sort.append('-n')

                            if self.get_option('compression-level'):
                                samtools_sort.extend(
                                    ['-l', self.get_option('compression-level')]
                                )

                            if self.get_option('max-mem-per-thread'):
                                samtools_sort.extend(
                                    ['-m', self.get_option('max-mem-per-thread')]
                                )

                            samtools_sort.extend(
                                ['-T', os.path.join(
                                    self.get_option('temp-sort-dir'),
                                     run_id), 
                                 '-',
                                 '-@', str(self.get_cores())]
                            )

                            pipe.add_command(samtools_sort)
                           
                            # 4. command: dd
                            dd_out = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize')
                            ]
                            pipe.add_command(
                                dd_out,
                                stdout_path = run.add_output_file(
                                    'alignments',
                                    '%s.sorted.%s' % (run_id, self.get_option('output-format')),
                                    input_paths)
                            )
