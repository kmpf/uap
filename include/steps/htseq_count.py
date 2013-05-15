import sys
from abstract_step import *
import pipeline
import re
import process_pool
import yaml

class HtSeqCount(AbstractStep):
    
    def __init__(self, pipeline):
        super(HtSeqCount, self).__init__(pipeline)
        
        self.set_cores(2)
        
        self.add_connection('in/alignments', constraints = {'min_files_per_run': 1, 'max_files_per_run': 1})
        self.add_connection('in/features', constraints = {'total_files': 1} )
        self.add_connection('out/counts')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('htseq-count')
        self.require_tool('grep')
        self.require_tool('invertGood')
        self.require_tool('samtools')

    def setup_runs(self, complete_input_run_info, connection_info):
        
        if not 'mode' in self.options:
            self.options['mode'] = 'union'
        if not 'stranded' in self.options:
            self.options['stranded'] = 'yes'
        if not 'type' in self.options:
            self.options['type'] = 'exon'
        if not 'idattr' in self.options:
            self.options['idattr'] = 'gene_id'
        if not 'fix_segemehl_copd' in self.options:
            self.options['fix_segemehl_copd'] = False
            
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
        
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:

                cat4m = [self.tool('cat4m'), run_info['info']['alignments_path']]
                pigz = [self.tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                grep = [self.tool('grep'), '-v', "\t\\*\t"]
                invertGood = [self.tool('invertGood')]
                samtools = [self.tool('samtools'), 'view', '-h', '-']
                htseq_count = [self.tool('htseq-count')]
                for key in ('mode', 'stranded', 'type', 'idattr'):
                    htseq_count.extend(['--%s' % key, self.options[key]])
                htseq_count.extend(['-', run_info['info']['features_path']])
        
                pipeline.append(cat4m)
                if run_info['info']['alignments_path'][-7:] == '.sam.gz':
                    pipeline.append(pigz)
                elif run_info['info']['alignments_path'][-4:] == '.bam':
                    pipeline.append(samtools)
                if self.options['fix_segemehl_copd'] == True:
                    pipeline.append(grep)
                    pipeline.append(invertGood)
                pipeline.append(htseq_count, stdout_path = run_info['info']['counts_path'])
        