import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml


class CuffLinks(AbstractStep):
    
    def __init__(self, pipeline):
        super(CuffLinks, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/features')
        self.add_connection('out/fpkm_tracking')
        self.add_connection('out/log')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('cufflinks')
        
    def setup_runs(self, complete_input_run_info, connection_info):
        if not 'library_type' in self.options:
            raise StandardError("Required option is missing: library_type.")
        
        output_run_info = dict()
            
        for run_id, info in connection_info['in/alignments']['runs'].items():
            transcripts_path = '%s-transcripts.gtf' % run_id
            skipped_path = '%s-skipped.gtf' % run_id
            genes_fpkm_tracking_path = '%s-genes.fpkm_tracking' % run_id
            isoforms_fpkm_tracking_path = '%s-isoforms.fpkm_tracking' % run_id
            log_path = '%s-cufflinks-log.txt' % run_id
            alignments_path = info.values()[0][0]
            run_info = {
                'output_files': {
                    'features': {
                        transcripts_path: [alignments_path],
                        skipped_path: [alignments_path]
                    },
                    'fpkm_tracking': {
                        genes_fpkm_tracking_path: [alignments_path],
                        isoforms_fpkm_tracking_path: [alignments_path]
                    },
                    'log': {
                        log_path: [alignments_path]
                    },
                },
                'info': {
                    'in-alignments': alignments_path,
                    'out-transcripts': transcripts_path,
                    'out-skipped': skipped_path,
                    'out-genes_fpkm_tracking': genes_fpkm_tracking_path,
                    'out-isoforms_fpkm_tracking': isoforms_fpkm_tracking_path
                }
            }
            output_run_info[run_id] = run_info
            
        return output_run_info

    def execute(self, run_id, run_info):
        
        cufflinks_out_path = self.get_temporary_path('cufflinks-out')
        
        with process_pool.ProcessPool(self) as pool:
            
            cufflinks = [
                self.tool('cufflinks'),
                '-o', cufflinks_out_path,
                '-p', '6', 
                "--library-type=%s" % self.options['library_type'],
                run_info['info']['in-alignments']
            ]

            pool.launch(cufflinks, stderr_path = run_info['output_files']['log'].keys()[0], 
            hints = {'writes': [
                run_info['output_files']['features'].keys()[0], 
                run_info['output_files']['features'].keys()[1], 
                run_info['output_files']['fpkm_tracking'].keys()[0],
                run_info['output_files']['fpkm_tracking'].keys()[0]
            ]})
            
        os.rename(os.path.join(cufflinks_out_path, 'transcripts.gtf'), run_info['info']['out-transcripts'])
        os.rename(os.path.join(cufflinks_out_path, 'skipped.gtf'), run_info['info']['out-skipped'])
        os.rename(os.path.join(cufflinks_out_path, 'genes.fpkm_tracking'), run_info['info']['out-genes_fpkm_tracking'])
        os.rename(os.path.join(cufflinks_out_path, 'isoforms.fpkm_tracking'), run_info['info']['out-isoforms_fpkm_tracking'])
        
        try:
            os.rmdir(cufflinks_out_path)
        except OSError:
            pass
