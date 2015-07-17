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
        self.require_tool('samtools')

        self.add_option('mode', str, default='mode')
        self.add_option('stranded', str, optional=False)
        self.add_option('type', str, default='exon')
        self.add_option('idattr', str, default='gene_id')
        
#        # TODO: remove fix_segemehl_copd option
#        self.add_option('fix_segemehl_copd', bool, default = False)
        
        
    def declare_runs(self):
        
        features_path = self.get_single_input_file_for_connection('in/features')

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            alignments_path = input_paths[0]

            with self.declare_run(run_id) as run:
                run.add_private_info('alignments_path', alignments_path)
                run.add_private_info('features_path', features_path)

                run.add_output_file('counts', '%s-counts.txt' % run_id, [alignments_path, features_path])
                run.new_exec_group()
#        features_path = connection_info['in/features']['runs'].values()[0].values()[0][0]
#        for run_id, info in connection_info['in/alignments']['runs'].items():
#            counts_path = '%s-counts.txt' % run_id
#            alignments_path = info.values()[0][0]
#            run_info = {
#                'output_files': {
#                    'counts': {
#                        counts_path: [alignments_path, features_path]
#                    }
#                },
#                'info': {
#                    'counts_path': counts_path,
#                    'alignments_path': alignments_path,
#                    'features_path': features_path
#                }
#            }
#            output_run_info[run_id] = run_info
        
#        return output_run_info
    
    
    def execute(self, run_id, run):
        
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                alignments_path = run.get_private_info('alignments_path')
                features_path = run.get_private_info('features_path')

                cat4m = [self.get_tool('cat4m'), alignments_path]
                pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1', '--stdout']
#                grep = [self.get_tool('grep'), '-v', "\t\\*\t"]
                samtools = [self.get_tool('samtools'), 'view', '-h', '-']
                htseq_count = [self.get_tool('htseq-count')]
                for key in ('mode', 'stranded', 'type', 'idattr'):
                    htseq_count.extend(['--%s' % key, self.get_option(key)])
                htseq_count.extend(['-', features_path])
        
                pipeline.append(cat4m)
                
                if alignments_path[-7:] == '.sam.gz':
                    pipeline.append(pigz)
                elif alignments_path[-4:] == '.bam':
                    pipeline.append(samtools)
                
                pipeline.append(htseq_count, stdout_path = run.get_single_output_file_for_annotation('counts'))
