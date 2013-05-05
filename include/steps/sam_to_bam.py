import sys
from abstract_step import *
import process_pool
import yaml


class SamToBam(AbstractStep):

    def __init__(self, pipeline):
        super(SamToBam, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('pigz')

    def setup_runs(self, complete_input_run_info, connection_info):
        # make sure files are available
        for key in ['genome']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find %s file: %s" % (key, self.options[key]))

        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for run_id, input_run_info in step_input_info.items():
                output_run_info[run_id] = {}
                output_run_info[run_id]['output_files'] = {}
                output_run_info[run_id]['output_files']['alignments']  = {}
                output_run_info[run_id]['output_files']['alignments'][run_id + '.bam'] = input_run_info['output_files']['alignments'].keys()
                output_run_info[run_id]['output_files']['alignments'][run_id + '.bam.bai'] = input_run_info['output_files']['alignments'].keys()
                output_run_info[run_id]['info'] = {}
                if len(input_run_info['output_files']['alignments'].keys()) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                output_run_info[run_id]['info']['in-sam']  = input_run_info['output_files']['alignments'].keys()[0]
                output_run_info[run_id]['info']['out-bam']  = run_id + '.bam'
                output_run_info[run_id]['info']['out-bai']  = run_id + '.bam.bai'

        return output_run_info

    def execute(self, run_id, run_info):
        sam_path = run_info['info']['in-sam']
        sorted_bam_path = run_info['info']['out-bam']
        sorted_bai_path = run_info['info']['out-bai']
        unsorted_bam_path = self.get_temporary_path('sam_to_bam_unsorted', 'output')
        
        # samtools view

        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), sam_path]
                pigz1 = [self.tool('pigz'), '--processes', '2', '-d', '-c']
                samtools = [self.tool('samtools'), 'view', '-Sbt', self.options['genome'], '-']
                
                pipeline.append(cat4m)
                pipeline.append(pigz1)
                pipeline.append(samtools, stdout_path = unsorted_bam_path)
        
        
        # samtools sort

        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat4m = [self.tool('cat4m'), unsorted_bam_path]
                samtools = [self.tool('samtools'), 'sort', '-', sorted_bam_path[:-4]]
                
                pipeline.append(cat4m)
                pipeline.append(samtools, hints = {'writes': [sorted_bam_path]})

        # samtools index
        
        with process_pool.ProcessPool(self) as pool:
            pool.launch([self.tool('samtools'), 'index', sorted_bam_path, '/dev/stdout'],
                stdout_path = sorted_bai_path, hints = {'reads': [sorted_bam_path]})

        os.unlink(unsorted_bam_path)