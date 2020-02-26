from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SamToSortedBam(AbstractStep):
    '''
    The step sam_to_sorted_bam builds on 'samtools sort' to sort SAM files and
    output BAM files.

    Sort alignments by leftmost coordinates, or by read name when -n is used.
    An appropriate @HD-SO sort order header tag will be added or an existing
    one updated if necessary.

    Documentation::

        http://www.htslib.org/doc/samtools.html
    '''

    def __init__(self, pipeline):
        super(SamToSortedBam, self).__init__(pipeline)

        self.set_cores(8)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')

        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        # Step was tested for samtools release 1.3.1
        self.require_tool('samtools')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

        self.add_option('sort-by-name', bool, default = False)
        self.add_option('genome-faidx', str, optional = False)
        self.add_option('temp-sort-dir', str, optional = True,
                        description = 'Intermediate sort files are stored into'
                        'this directory.')

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "256k")

    def runs(self, run_ids_connections_files):

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
                    sortpath = os.path.abspath(self.get_option('temp-sort-dir'))
                    if not os.path.isdir(sortpath):
                        #dir not present
                        raise StepError(self, "Directory %s not found" % self.get_option('temp-sort-dir'))
                    if not os.access(sortpath, os.W_OK):
                        #not accessible
                        raise StepError(self, "Directory %s not accessible." % self.get_option('temp-sort-dir'))
                else:
                    sortpath =  (run.get_output_directory_du_jour_placeholder()  + '/')


                with run.new_exec_group() as exec_group:

                    with exec_group.add_pipeline() as pipe:
                        # 1. command: Read file in 4MB chunks
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

                        # 2. command: Convert sam to bam
                        samtools_view = [
                            self.get_tool('samtools'), 'view',
                            '-S', '-b', '-t',
                            self.get_option('genome-faidx'), '-',
                            '-@', '1'
                        ]
                        pipe.add_command(samtools_view)

                        # 3. command: Sort BAM input
                        samtools_sort = [
                            self.get_tool('samtools'), 'sort',
                            '-O', 'bam'
                        ]
                        if self.get_option('sort-by-name'):
                            samtools_sort.append('-n')
                        samtools_sort.extend(
                            ['-T',
                             os.path.join(sortpath, run_id),
                             '-',
                             '-@', '6']
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
                                '%s.sorted.bam' % run_id,
                                input_paths)
                        )
