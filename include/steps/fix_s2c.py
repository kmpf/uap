import sys
from abstract_step import *
import process_pool
import yaml


class S2CFix(AbstractStep):

    def __init__(self, pipeline):
        super(S2CFix, self).__init__(pipeline)

        
        self.set_cores(2)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')

        
        self.require_tool('fix_s2c')
        self.require_tool('pigz')
        self.require_tool('cat4m')
        self.require_tool('samtools')

    def setup_runs(self, complete_input_run_info, connection_info):
        output_run_info = {}
        
        for run_id, info in connection_info['in/alignments']['runs'].items():
            fix_path = '%s-s2c-fixed.bam' % run_id

            alignments_path = info.values()[0][0]
            run_info = {
                'output_files': {
                    'alignments': {
                        fix_path: [alignments_path]
                    }

                },
                'info': {
                    'fix_path': fix_path,
                    'alignments_path': alignments_path
                }
            }
            output_run_info[run_id] = run_info
        
        return output_run_info


    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), run_info['info']['alignments_path']]
                samtools = [self.tool('samtools'), 'view', '-h', '-']
                fix_s2c = [self.tool('fix_s2c'), ]
                samtools_2 = [self.tool('samtools'), 'view', '-Shb', '-']
                samtools_sort = [self.tool('samtools'), 'sort', '-', run_info['info']['fix_path'][:-4]]

                
                pipeline.append(cat4m)
                pipeline.append(samtools)
                pipeline.append(fix_s2c)
                pipeline.append(samtools_2)
                pipeline.append(samtools_sort, hints = {'writes': [run_info['info']['fix_path']]})
