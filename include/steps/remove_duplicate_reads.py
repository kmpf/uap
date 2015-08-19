import sys
import yaml

from ..abstract_step import *
from .. import misc
from .. import process_pool

class RemoveDuplicates(AbstractStep):


    def __init__(self, pipeline):
        super(RemoveDuplicates, self).__init__(pipeline)

        self.set_cores(12)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/metrics')
        
        self.require_tool('cat')
        self.require_tool('MarkDuplicates')

#        self.add_option('index', str)

    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                run.new_exec_group()
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")
                basename_parts = os.path.basename(input_paths[0]).split('.')
                if all(suffix not in ['sam', 'bam'] for suffix in basename_parts):
                    raise StandardError("The file %s seems not to be a SAM or BAM "
                                        "file. At least the suffix is wrong." % input_paths[0])

                run.add_private_info('in-alignment', input_paths[0])
                run.add_output_file('alignments', '%s-rm-dup.bam' % run_id,
                                    input_paths)
                run.add_output_file('metrics', '%s-rm-dup-metrics.txt' % run_id, input_paths)

    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            alignment_path = run.get_private_info('in-alignment')
            with pool.Pipeline(pool) as pipeline:
#                cat_input = [self.get_tool('cat'), alignment_path]
                mark_duplicates = [
                    self.get_tool('MarkDuplicates'),
                    'INPUT=%s' % alignment_path,
                    'OUTPUT=%s' % run.get_single_output_file_for_annotation('alignments'),
                    'METRICS_FILE=%s' % run.get_single_output_file_for_annotation('metrics'),
                    'REMOVE_DUPLICATES=true'                
                    ]
                                                  
#                cat_output = [self.get_tool('cat'), '-']
#                print(cat_input)
#                print(mark_duplicates)
#                print(cat_output)

#                pipeline.append(cat_input)
                pipeline.append(mark_duplicates, 
                                hints={'writes': 
                                       [run.get_single_output_file_for_annotation('metrics'),
                                        run.get_single_output_file_for_annotation('alignments')]
                                       })
#                pipeline.append(cat_output, 
#                                stdout_path = )

