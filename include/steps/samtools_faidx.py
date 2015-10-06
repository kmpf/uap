import os

from abstract_step import AbstractStep

class SamtoolsFaidx(AbstractStep):

    def __init__(self, pipeline):
        super(SamtoolsFaidx, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/sequence')
        self.add_connection('out/indices')
        
        self.require_tool('samtools')
        self.require_tool('mv')
        
    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/sequence']

                if input_paths == [None]:
                    run.add_empty_output_connection("sequence")
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one sequence file.")
                elif os.path.splitext(os.path.basename(input_paths[0]))[1] \
                not in ['.fa', '.fna', '.fasta']:
                    raise StandardError("The input file %s does not seem to be "
                                        "a FASTA file." % input_paths[0])
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