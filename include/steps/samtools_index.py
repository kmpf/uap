import sys
import os
from abstract_step import AbstractStep
import yaml

class SamtoolsIndex(AbstractStep):

    def __init__(self, pipeline):
        super(SamToBam, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/indices')
        self.add_connection('out/index_stats')
        
        self.require_tool('ln')
        self.require_tool('samtools')

        self.add_option('index_type', str, choices = ['bai', 'csi'],
                        optional = False)

    def runs(self, run_ids_connections_files):
                
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                # Add empty out connection if we have an empty in connection
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                    run.add_empty_output_connection("indices")
                # Fail if we haven't exactly one input file
                elif len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                # Fail if the input is not a bam file
                elif os.path.splitext(input_paths[0])[1] not in ['.bam']:
                    raise StandardError(
                        "The file %s seems not to be a BAM file. At "
                        "least the suffix is wrong." % input_paths[0]
                    )
                # Everything seems fine, lets start
                else:
                    input_bam = input_paths[0]
                    # At first create the index and a symlink to original BAM
                    with run.new_exec_group() as index_exgr:
                        # 1. command: Create symbolic link to original bam file
                        # (use absolute path)
                        ln = [self.get_tool('ln'), '-s', input_bam]
                        index_exgr.add_command(ln)
                        # 2. command: Index bam file
                        samtools_index = [self.get_tool('samtools'), 'index']
                        if self.get_option('index_type') == 'bai':
                            samtools_index.append('-b')
                            run.add_output_file(
                                'indices',
                                '%s.bai' % os.path.basename(input_bam),
                                input_paths
                            )
                        elif self.get_option('index_tpye') == 'csi':
                            samtools_index.append('-c')
                            run.add_output_file(
                                'indices',
                                '%s.csi' % os.path.basename(input_bam),
                                input_paths
                            )
                        samtools_index.append(input_bam)
                        index_exgr.add_command(samtools_index)
                    # Calculate samtools idxstats
                    with run.new_exec_group() as idxstats_exgr:
                        samtools_idxstats = [
                            self.get_tool('samtools'), 'idxstats']
                        samtools_idxstats.extend(input_bam)
                        idxstats_exgr.add_command(
                            samtools_idxstats,
                            stdout_path = run.add_output_file(
                                'index_stats',
                                '%s_idxstats.txt' % os.path.basename(input_bam),
                                input_paths
                            )
                        )
