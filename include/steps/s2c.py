import sys
from abstract_step import *
import unix_pipeline
import yaml


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
        self.require_tool('cat4m')
        
    def setup_runs(self, complete_input_run_info, connection_info):
        output_run_info = {}
        
        for run_id, info in connection_info['in/alignments']['runs'].items():
            s2c_path = '%s-cufflinks-compatible-sorted.bam' % run_id
            log_path = '%s-log.txt' % run_id
            alignments_path = info.values()[0][0]
            run_info = {
                'output_files': {
                    'alignments': {
                        s2c_path: [alignments_path]
                    },
                    'log': {
                        log_path: [alignments_path]
                    },
                },
                'info': {
                    's2c_path': s2c_path,
                    'log_path': log_path,
                    'alignments_path': alignments_path
                }
            }
            output_run_info[run_id] = run_info
        
        return output_run_info

    def execute(self, run_id, run_info):
        p = unix_pipeline.UnixPipeline()
        
        cat4m = [self.tool('cat4m'), run_info['info']['alignments_path']]
        pigz = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
        s2c = [self.tool('s2c'), '-s', '/dev/stdin', '-o', self._temp_directory]
        samtools = [self.tool('samtools'), 'view', '-Sb', '-']
        samtools_sort = [self.tool('samtools'), 'sort', '-', run_info['info']['s2c_path'][:-4]]
        
        p.append(cat4m)
        p.append(pigz)
        p.append(s2c, stderr_path = run_info['info']['log_path'])
        p.append(samtools)
        p.append(samtools_sort)
                
        unix_pipeline.wait()

        p = unix_pipeline.UnixPipeline()
        