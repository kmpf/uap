import sys
import yaml

from ..abstract_step import *
from .. import process_pool

class S2C(AbstractStep):

    def __init__(self, pipeline):
        super(S2C, self).__init__(pipeline)
        
        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log')
        
        self.require_tool('s2c')
        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cat')
        
    def declare_runs(self):

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file., but got this %s" % input_paths)
                run.add_private_info('in-alignment', input_paths[0])
                run.add_output_file('alignments', '%s-cufflinks-compatible-sorted.bam' % run_id, input_paths)
                run.add_output_file('log', '%s-log.txt' % run_id, input_paths)
                run.new_exec_group()
#        for run_id, info in connection_info['in/alignments']['runs'].items():
#            s2c_path = '%s-cufflinks-compatible-sorted.bam' % run_id
#            log_path = '%s-log.txt' % run_id
#            alignments_path = info.values()[0][0]
#            run_info = {
#                'output_files': {
#                    'alignments': {
#                        s2c_path: [alignments_path]
#                    },
#                    'log': {
#                        log_path: [alignments_path]
#                    },
#                },
#                'info': {
#                    's2c_path': s2c_path,
#                    'log_path': log_path,
#                    'alignments_path': alignments_path
#                }
#            }
#            output_run_info[run_id] = run_info
#        
#        return output_run_info

    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                alignments_path = run.get_private_info('in-alignment')
                cat = [self.get_tool('cat'), alignments_path]
                pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                s2c = [self.get_tool('s2c'), '-s', '/dev/stdin', '-o', self._temp_directory]
                samtools = [self.get_tool('samtools'), 'view', '-Sb', '-']
                samtools_sort = [self.get_tool('samtools'), 'sort', '-', run.get_single_output_file_for_annotation('alignments')[:-4]]
                
                pipeline.append(cat)
                pipeline.append(pigz)
                pipeline.append(s2c, stderr_path = run.get_single_output_file_for_annotation('log'))
                pipeline.append(samtools)
                pipeline.append(samtools_sort, hints = {'writes': run.get_single_output_file_for_annotation('alignments')})
