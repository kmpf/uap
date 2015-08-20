import sys

from ..abstract_step import *
from .. import process_pool
from .. import yaml

class S2cFix(AbstractStep):

    def __init__(self, pipeline):
        super(S2cFix, self).__init__(pipeline)
        
        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')

        
        self.require_tool('fix_s2c')
        self.require_tool('pigz')
        self.require_tool('cat')
        self.require_tool('samtools')

    def declare_runs(self):

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                if len(input_paths) != 1:
                    raise Exception("Expected exactly one alignments file, but got this %s" % input_paths)
                run.add_output_file('alignments', '%s-s2c-fixed.bam' % run_id, input_paths)
                run.add_private_info('in-alignments', input_paths[0])
                run.new_exec_group()

#        for run_id, info in connection_info['in/alignments']['runs'].items():
#            fix_path = 

#            alignments_path = info.values()[0][0]
#            run_info = {
#                'output_files': {
#                    'alignments': {
#                        fix_path: [alignments_path]
#                    }
#                },
#                'info': {
#                    'fix_path': fix_path,
#                    'alignments_path': alignments_path
#                }
#            }
#            output_run_info[run_id] = run_info
#        return output_run_info

    def execute(self, run_id, run_info):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                in_alignment = run.get_private_info('in-alignments')
                cat = [self.get_tool('cat'), in_alignment]
                samtools = [self.get_tool('samtools'), 'view', '-h', '-']
                fix_s2c = [self.get_tool('fix_s2c'), ]
                samtools_2 = [self.get_tool('samtools'), 'view', '-Shb', '-']
                samtools_sort = [self.get_tool('samtools'), 'sort','-n', '-', run.get_single_output_file_for_annotation('alignments')[:-4]]

                
                pipeline.append(cat)
                pipeline.append(samtools)
                pipeline.append(fix_s2c)
                pipeline.append(samtools_2)

                pipeline.append(samtools_sort, hints = {'writes': run.get_single_output_file_for_annotation('alignments')})
