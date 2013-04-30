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
        if connection_info['in/features']['counts']['total_files'] != 1:
            raise StandardError("Expected a single in/features file:\n%s" % yaml.dump(connection_info, default_flow_style = False))
        
        output_run_info = {}
        
        for run_id, info in connection_info['in/alignments']['runs'].items():
            counts_path = '%s-counts.txt' % run_id
            if len(info.values()) != 1:
                raise StandardError("")
            if len(info.values()[0]) != 1:
                raise StandardError("")
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
        return
        # basic sanity check
        if len(run_info['output_files']['reads']) != 1:
            raise StandardError("Expected a single output file.")

        # set up processes
        cat4m = [self.tool('cat4m')]
        cat4m.extend(*sorted(run_info['output_files']['reads'].values()))

        pigz1 = [self.tool('pigz'), '--processes', '1', '--decompress', '--stdout']

        cutadapt = [self.tool('cutadapt'), '-a', run_info['info']['adapter'], '-']

        pigz2 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '3', '--stdout']

        # create the pipeline and run it
        p = unix_pipeline.UnixPipeline()
        p.append(cat4m)
        p.append(pigz1)
        p.append(cutadapt, stderr_path = run_info['output_files']['log'].keys()[0])
        p.append(pigz2, stdout_path = run_info['output_files']['reads'].keys()[0])

        p.seal()
        
        unix_pipeline.wait()
