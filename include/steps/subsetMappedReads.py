from uaperrors import UAPError
import sys
from abstract_step import *
import process_pool
import yaml
import os
from logging import getLogger

logger=getLogger('uap_logger')

class subsetMappedReads(AbstractStep):

    '''
    subsetMappedReads selects a provided number of mapped reads from a
    file in .sam or .bam format. Depending on the set options the
    first N mapped reads and their mates (for paired end sequencing)
    are returned in .sam format. If the number of requested reads
    exceeds the number of available mapped reads, all mapped reads are
    returned. 
    '''

    def __init__(self, pipeline):
        super(subsetMappedReads, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log')

        self.require_tool('samtools')
        self.require_tool('dd')
        self.require_tool('pigz')
        self.require_tool('head')
        self.require_tool('cat')
        
        self.add_option('genome-faidx', str, optional = False)
        self.add_option('Nreads', str, optional=False,
                        description='Number of reads to extract from input file. ')
        self.add_option('paired_end', bool, optional = False,
                        description='The reads are expected to have a mate, due to paired end sequencing.')

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

                if not self.is_option_set_in_config('Nreads'):
                    raise UAPError("Required option 'Nreads' not set in your configuration file")

                with run.new_exec_group() as exec_group:

                    with exec_group.add_pipeline() as pipe:
                        # 1. command: Read file in 4MB chunks
                        dd_in = [self.get_tool('dd'),
                                 'ibs=4M',
                                 'if=%s' % input_paths[0]]
                        pipe.add_command(dd_in)

                        # 1.1 command: Uncompress file to fifo
                        if is_gzipped:
                            pigz = [self.get_tool('pigz'),
                                    '--decompress',
                                    '--processes', str(self.get_cores()),
                                    '--stdout'
                            ]
                            pipe.add_command(pigz)

                        # 2. command: Read sam file
                        # extract only reads that were aligned and include only pairs
                        if self.get_option('paired_end'):
                            samtools_view = [
                                self.get_tool('samtools'), 'view', '-F', '0x04', '-f', '0x03', '-h',
#                                '-t', os.path.abspath(self.get_option('genome-faidx')),
                                '-'
                            ]
                        else:
                            samtools_view = [
                                self.get_tool('samtools'), 'view', '-F', '0x04', '-h',
#                                self.get_tool('samtools'), 'view', '-c',
#                                '-t', os.path.abspath(self.get_option('genome-faidx')),
                                '-'
                            ]
                        pipe.add_command(samtools_view)

                        # 3. extract the first Nreads
                        # include mate and the lines for the header of the sam format
                        N = int(self.get_option('Nreads')) * 2 + 24
                        get_Nreads = [
                            self.get_tool('head'), '-%s' % N
                        ]
                        pipe.add_command(get_Nreads)

                        dd_out = [self.get_tool('dd'), 'obs=4M']
                        pipe.add_command(
                            dd_out,
                            stdout_path = run.add_output_file(
                                'alignments',
                                '%s.N%s.reads.sam' % (run_id, self.get_option('Nreads')),
                                input_paths)
                        )
