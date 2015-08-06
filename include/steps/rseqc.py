import sys
from abstract_step import *
import process_pool
import yaml

class RSeQC(AbstractStep):
    '''
    The RSeQC step can be used to evaluate aligned reads in a BAM file. RSeQC does not only report raw sequence-based metrics, but also quality control metrics like read distribution, gene coverage, and sequencing depth.
    '''

    def __init__(self, pipeline):
        super(RSeQC, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/alignments')
        self.add_connection('out/bam_stat')
#	self.add_connection('out/junction_annotation')
	self.add_connection('out/infer_experiment')
	self.add_connection('out/read_distribution')
#	self.add_connection('out/junction_saturation')
#	self.add_connection('out/RPKM_saturation')
                
        self.require_tool('cat')
        self.require_tool('bam_stat.py')
#	self.require_tool('junction_annotation.py')
	self.require_tool('infer_experiment.py')
	self.require_tool('read_distribution.py')
#	self.require_tool('junction_saturation.py')
#	self.require_tool('RPKM_saturation.py')

	self.add_option('reference', str)
        
    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:

                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignment file.")

                basename = os.path.basename(input_paths[0]).split('.')[0]
                run.add_output_file('bam_stat', basename + '.bam_stats.txt', input_paths)
		run.add_output_file('infer_experiment', basename + '.infer_experiment.txt', input_paths)
		run.add_output_file('read_distribution', basename + '.read_distribution.txt', input_paths)
                run.add_private_info('in-bam', input_paths[0])

    def execute(self, run_id, run):
        bam_in_path = run.get_private_info('in-bam')
        out_path_bam_stat = run.get_single_output_file_for_annotation('bam_stat')
	out_path_infer_experiment = run.get_single_output_file_for_annotation('infer_experiment')
	out_path_read_distribution = run.get_single_output_file_for_annotation('read_distribution')
	ref = self.get_option('reference')
	
        with process_pool.ProcessPool(self) as pool:
		pool.launch([self.get_tool('bam_stat.py'), '-i', bam_in_path], stderr_path = out_path_bam_stat)  # bam_stat.py writes the output to stderr
		
	with process_pool.ProcessPool(self) as pool:
		pool.launch([self.get_tool('infer_experiment.py'), '-i', bam_in_path, '-r', ref], stdout_path = out_path_infer_experiment)
		
	with process_pool.ProcessPool(self) as pool:
		pool.launch([self.get_tool('read_distribution.py'), '-i', bam_in_path, '-r', ref], stdout_path = out_path_read_distribution)





