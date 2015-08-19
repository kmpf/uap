import sys
import re
import yaml

from ..abstract_step import *
from .. import process_pool

class BedgraphToBigwig(AbstractStep):

    def __init__(self, pipeline):
        super(BedgraphToBigwig, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/tracks')
        self.add_connection('out/tracks')
        
        self.require_tool('bedGraphToBigWig')

        self.add_option('genome', str, optional=False)

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/tracks'):
            with self.declare_run(run_id) as run:
                bigwig_files = list()
                for input_file in input_paths:
                    validated_input_files = list()
                    for suffix in ['.bg', '.bedgraph']:
                        if input_file[-len(suffix):] == suffix:
                            validated_input_files.append(input_file)
                    if len(validated_input_files) == 1:
                        bigwig_file = os.path.basename(input_file)[:-len(suffix)] + '.bw'
                        run.add_output_file('tracks', bigwig_file, validated_input_files)
                        bigwig_files.append(bigwig_file)
                    else:
                        raise StandardError("%s file suffix is not '%s'. Please provide a BEDGRAPH file" % (input_file, suffix))
                run.add_private_info('output_files', bigwig_files)
                run.new_exec_group()
                    
    def execute(self, run_id, run):
        output_files = run.get_output_files_abspath()
        for bigwig_file, input_files in output_files['tracks'].items():
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    bedgraph_to_bigwig = [
                        self.get_tool('bedGraphToBigWig'),
                        input_files[0],
                        self.get_option('genome'),
                        bigwig_file]
                    pipeline.append(bedgraph_to_bigwig)

