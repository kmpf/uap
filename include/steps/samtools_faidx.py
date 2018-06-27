import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SamtoolsFaidx(AbstractStep):
    '''
    Index reference sequence in the FASTA format or extract subsequence from
    indexed reference sequence. If no region is specified, faidx will index the
    file and create <ref.fasta>.fai on the disk. If regions are specified, the
    subsequences will be retrieved and printed to stdout in the FASTA format.

    The sequences in the input file should all have different names. If they do
    not, indexing will emit a warning about duplicate sequences and retrieval
    will only produce subsequences from the first sequence with the duplicated
    name. 
    '''

    def __init__(self, pipeline):
        super(SamtoolsFaidx, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/sequence')
        self.add_connection('out/indices')
        
        self.require_tool('samtools')
        self.require_tool('mv')

        # Use of optional '[region1 [...]]' argument is not supported, as it
        # changes the behaviour of the command
        
    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/sequence']

                if input_paths == [None]:
                    run.add_empty_output_connection("sequence")
                elif len(input_paths) != 1:
                    logger.error("Expected exactly one sequence file.")
                    sys.exit(1)
                elif os.path.splitext(os.path.basename(input_paths[0]))[1] \
                not in ['.fa', '.fna', '.fasta']:
                    logger.error("The input file %s does not seem to be "
                                 "a FASTA file." % input_paths[0])
                    sys.exit(1)
                else:
                    with run.new_exec_group() as faidx_group:
                        samtools_faidx = [
                            self.get_tool('samtools'),
                            'faidx', input_paths[0]]
                        faidx_group.add_command(samtools_faidx)

                    with run.new_exec_group() as mv_group:
                        mv = [self.get_tool('mv'),
                              input_paths[0] + '.fai',
                              run.add_output_file(
                                  "indices",
                                  "%s.fai" % os.path.basename(input_paths[0]),
                                  input_paths
                              )
                        ]
                        mv_group.add_command(mv)
