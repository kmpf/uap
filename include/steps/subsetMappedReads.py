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
    first N mapped reads or a random selection (not implemented yet)
    of N mapped reads are returned in .sam format. If the number of
    requested reads exceeds the number of available mapped reads, all
    mapped reads are returned.  NOTE: the random version is not yet
    implemented - don't use it!
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
        self.require_tool('head') # do we really need this to require?

        self.add_option('genome-faidx', str, optional = False)
        self.add_option('Nreads', str, optional=False,
                        description='Number of reads to extract from input file.')
        # any idea to do this fast?
        # - shuffle the number of ALL mapped reads, return the first N
        # - generate N (distinct!) random numbers, return the appropriate reads
        #self.add_option('random', bool, default=False, optional=True,
                        #description='Select the reads randomly.')

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    logger.error("Expected exactly one alignments file.")
                    sys.exit(1)
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                if not self.is_option_set_in_config('Nreads'):
                    logger.error("Required option 'Nreads' not set in your configuration file")
                    sys.exit(1)

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
                                    '--processes', '1',
                                    '--stdout']
                            pipe.add_command(pigz)

                        # 2. command: Read sam file
                        samtools_view = [
                            self.get_tool('samtools'), 'view',
                            '-S', '-t', self.get_option('genome-faidx'),
                            '-'
                        ]
                        pipe.add_command(samtools_view)

                        # 3. extract the first Nreads
                        get_Nreads = [
                            self.get_tool('head'), '-%s' % self.get_option('Nreads')
                        ]
                        pipe.add_command(get_Nreads)

                        # 4. command: Write sam file
                        samtools_write = [
                            self.get_tool('samtools'), 'view',
                            '-S', '-'
                        ]
                        pipe.add_command(samtools_write)

                        # 5. command: dd
                        dd_out = [self.get_tool('dd'), 'obs=4M']
                        pipe.add_command(
                            dd_out,
                            stdout_path = run.add_output_file(
                                'alignments',
                                '%s.N%s.reads.sam' % (run_id, self.get_option('Nreads')),
                                input_paths)
                        )
