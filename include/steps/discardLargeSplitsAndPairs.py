import sys
from abstract_step import *
import process_pool
import yaml
import os
from logging import getLogger

logger=getLogger('uap_logger')

class discardLargeSplitsAndPairs (AbstractStep):

    '''
    discardLargeSplitsAndPairs reads SAM formatted alignments of the
    mapped reads. It discards all split reads that skip more than
    splits_N nucleotides in their alignment to the ref genome. In
    addition, all read pairs that are mapped to distant region such
    that the final template will exceed N_mates nucleotides will also
    be discarded. All remaining reads are returned in SAM format. The
    discarded reads are also collected in a SAM formatted file and a
    statistic is returned.
    '''

    def __init__(self, pipeline):
        super(discardLargeSplitsAndPairs, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/alignments')  # sam aln
        self.add_connection('out/alignments') # contains the kept reads
        self.add_connection('out/log')        # contains the discarded reads
        self.add_connection('out/stats')      # contains a statistic

        self.require_tool('dd')
        self.require_tool('pigz')
        self.require_tool('discardLargeSplitsAndPairs')
        
        self.add_option('N_splits', str, optional = False,
                        description='Size of the skipped region within a split read (in nucleotides). Split Reads that skip more nt than this value are discarded.')
        self.add_option('M_mates', str, optional=False,
                        description='Size of template (in nucleotides) that would arise from a read pair. Read pairs that exceed this value are discarded. ')

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

                        # 2. command: Process sam file
                        # create the names of the out connections
                        outfile = run.add_output_file('alignments',
                                                      '%s.reduced.sam' % run_id,
                                                      input_paths)
                        logfile = run.add_output_file('log',
                                                      '%s.discarded.sam' % run_id,
                                                      input_paths)
                        statsfile = run.add_output_file('stats',
                                                        '%s.statistics.txt' % run_id,
                                                        input_paths)
                        # construct cmd
                        discard_cmd = [self.get_tool('discardLargeSplitsAndPairs'),
                                       '--N_splits', self.get_option('N_splits'),
                                       '--M_mates', self.get_option('M_mates'),
                                       '--logfile', logfile,
                                       '--statsfile', statsfile
                                       ]
                        # execute cmd                              
                        pipe.add_command(discard_cmd, stdout_path = outfile)
