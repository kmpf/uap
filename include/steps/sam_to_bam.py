import sys
from abstract_step import *
import unix_pipeline
import yaml


class SamToBam(AbstractStep):
    def __init__(self, pipeline):
        super(SamToBam, self).__init__(pipeline)
        self.set_cores(4)

    def setup_runs(self, complete_input_run_info):
        # make sure tools are available
        self.tool('samtools')
        self.tool('pigz')
        self.tool('cat4m')

        # make sure files are available
        for key in ['genome']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find " + key + " file: " + self.options[key])

        output_run_info = {}
        for run_id, input_run_info in complete_input_run_info.items():
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
        unsorted_bam_path = self.get_temporary_path('sam_to_bam-unsorted-', '.bam')
        sorted_bam_path = run_info['info']['out-bam']
        sorted_bai_path = run_info['info']['out-bai']
        
        # samtools view

        cat4m = [self.tool('cat4m'), sam_path]
        pigz1 = [self.tool('pigz'), '--processes', '2', '-d', '-c']
        samtools = [self.tool('samtools'), 'view', '-Sbt', self.options['genome'], '-']
        
        p = unix_pipeline.create_pipeline()
        p.append(cat4m)
        p.append(pigz1)
        p.append(samtools, stdout = open(unsorted_bam_path, 'w'))
        
        unix_pipeline.wait()
        
        # samtools sort

        cat4m = [self.tool('cat4m'), unsorted_bam_path]
        samtools = [self.tool('samtools'), 'sort', '-', sorted_bam_path[:-4]]
        
        p = unix_pipeline.create_pipeline()
        p.append(cat4m)
        p.append(samtools)
        
        unix_pipeline.wait()

        # samtools index
        
        unix_pipeline.launch([self.tool('samtools'), 'index', sorted_bam_path, '/dev/stdout'],
            stdout = open(sorted_bai_path, 'w'))

        unix_pipeline.wait()
        
        os.unlink(unsorted_bam_path)