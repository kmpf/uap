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
        self.require_tool('grep')
        self.require_tool('invertGood')

    def setup_runs(self, complete_input_run_info, connection_info):
        
        output_run_info = {}
        
        features_path = connection_info['in/features']['runs'].values()[0].values()[0][0]
        for run_id, info in connection_info['in/alignments']['runs'].items():
            counts_path = '%s-counts.txt' % run_id
            alignments_path = info.values()[0][0]
            run_info = {
                'output_files': {
                    'counts': {
                        counts_path: [alignments_path, features_path]
                    }
                },
                'info': {
                    'counts_path': counts_path,
                    'alignments_path': alignments_path,
                    'features_path': features_path
                }
            }
            output_run_info[run_id] = run_info
        
        return output_run_info
    
    
    def execute(self, run_id, run_info):
        p = unix_pipeline.UnixPipeline()
        
        if not 'mode' in self.options:
            self.options['mode'] = 'union'
        if not 'stranded' in self.options:
            self.options['stranded'] = 'yes'
        if not 'type' in self.options:
            self.options['type'] = 'exon'
        if not 'idattr' in self.options:
            self.options['idattr'] = 'gene_id'

        cat4m = [self.tool('cat4m'), run_info['info']['alignments_path']]
        pigz = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
        grep = [self.tool('grep'), '-v', "\t\\*\t"]
        invertGood = [self.tool('invertGood')]
        htseq_count = [self.tool('htseq-count')]
        for key in ('mode', 'stranded', 'type', 'idattr'):
            htseq_count.extend(['--%s' % key, self.options[key]])
        htseq_count.extend(['-', run_info['info']['features_path']])
        
        p.append(cat4m)
        p.append(pigz)
        p.append(grep)
        p.append(invertGood)
        p.append(htseq_count, stdout_path = run_info['info']['counts_path'])
                
        unix_pipeline.wait()
        