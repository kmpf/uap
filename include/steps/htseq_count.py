import sys
from abstract_step import *
import pipeline
import re
import unix_pipeline
import yaml

class HtSeqCount(AbstractStep):
    
    def __init__(self, pipeline):
        super(HtSeqCount, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments', constraints = {'min_files_per_run': 1, 'max_files_per_run': 1})
        self.add_connection('in/features', constraints = {'total_files': 1} )
        self.add_connection('out/counts')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('htseq-count')

    def setup_runs(self, complete_input_run_info, connection_info):
        
        output_run_info = {}
        
        for run_id, info in connection_info['in/alignments']['runs'].items():
            counts_path = '%s-counts.txt' % run_id
            alignments_path = info.values()[0][0]
            run_info = {
                'output_files': {
                    'counts': {
                        counts_path: [alignments_path]
                    }
                },
                'info': {
                    'counts_path': counts_path,
                    'alignments_path': alignments_path,
                    'features_path': connection_info['in/features']['runs'].values()[0].values()[0][0]
                }
            }
            output_run_info[run_id] = run_info
        
        return output_run_info
    
    
    def execute(self, run_id, run_info):
        run_info['info']['features_path'] = '/home/michael/programming/rnaseq-pipeline/out/gencode-7898/genes.gtf.gz'

        features_fifo = unix_pipeline.mkfifo()
        
        cat4m = [self.tool('cat4m'), run_info['info']['features_path'], '-o', features_fifo]
        unix_pipeline.launch(cat4m)
        
        p = unix_pipeline.UnixPipeline()
        
        cat4m2 = [self.tool('cat4m'), run_info['info']['alignments_path']]
        htseq_count = [self.tool('htseq-count'), '-', features_fifo]
        
        p.append(cat4m2)
        p.append(htseq_count, stdout_path = run_info['info']['counts_path'])
                
        unix_pipeline.wait()
        
        os.unlink(features_fifo)
