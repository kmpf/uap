import sys
import os
from abstract_step import AbstractStep

class SamToSortedBam(AbstractStep):

    def __init__(self, pipeline):
        super(SamToSortedBam, self).__init__(pipeline)
        
        self.set_cores(8)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('pigz')

        self.add_option('sort-by-name', bool, default = False)
        self.add_option('genome-faidx', str, optional = False)
        self.add_option('temp-sort-directory', str, optional = False,
                        description = 'Intermediate sort files are stored into'
                        'this directory.')

    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                    with run.new_exec_group() as exec_group:

                        with exec_group.add_pipeline() as pipe:
                            # 1. command: Read file in 4MB chunks
                            dd_in = [self.get_tool('dd'),
                                     'ibs=4M',
                                     'if=%s' % input_paths[0]]
                            pipe.add_command(dd_in)

                            if is_gzipped:
                                # 1.1 command: Uncompress file to fifo
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
                                ['-T', os.path.join(
                                    self.get_option('temp-sort-directory'), run_id), 
                                 '-',
                                 '-@', '6']
                            )
                            pipe.add_command(samtools_sort)

                            # 4. command:
                            dd_out = [self.get_tool('dd'), 'obs=4M']
                            pipe.add_command(
                                dd_out,
                                stdout_path = run.add_output_file(
                                    'alignments',
                                    '%s.sorted.bam' % run_id,
                                    input_paths)
                            )
