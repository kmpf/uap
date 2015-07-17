import sys
from abstract_step import *
import process_pool
import yaml


class CountMappedReads(AbstractStep):

    def __init__(self, pipeline):
        super(CountMappedReads, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignment-counts')
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('pigz')

        self.add_option('set_FLAG_bits', list, default=0)
        self.add_option('unset_FLAG_bits', list, default=0)
        self.add_option('exclude_MAPQ_smaller_than', int, default=0)

    def declare_runs(self):
        # set_FLAG_bits and unset_FLAG_bits must have equal length,
        # because both are set everytime samtools is called
        set_bits = self.get_option('set_FLAG_bits')
        unset_bits = self.get_option('unset_FLAG_bits')
        
        if len(set_bits) != len(unset_bits):
            raise StandardError("set_FLAG_bits has %s elements but "
                                "unset_FLAG_bits has %s elements. Both "
                                "options need to have same length.")
        
        # Get run_id and input files for 'in/alignment' connection
        sample_alignment_paths_dict = dict()
        alignment_files = list()

        for run_id, alignment_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            # Only one BAM/SAM file expected as input file
            if len(alignment_paths) != 1:
                raise StandardError("Expected exactly one file with mapped reads.")
                
            sample_alignment_paths_dict[run_id] = alignment_paths[0]
            alignment_files.extend(alignment_paths)
                                
        # Name the run as the step
        run_id = self.get_step_name()

        with self.declare_run(run_id) as run:
            run.add_private_info('sample-alignment-paths', 
                                 sample_alignment_paths_dict)
            run.add_output_file('alignment-counts', 
                                '%s_alignment_counts.txt' % run_id,
                                alignment_files)
            run.new_exec_group()

    def execute(self, run_id, run):
        # Get the names of the announced output files
        out_alignment_counts_file = run.get_single_output_file_for_annotation('alignment-counts')

        # Get the input files per sample and connection
        sample_alignment_paths_dict = run.get_private_info('sample-alignment-paths')

        # Get names of temporary output directories
        temp_alignments_dir = self.get_temporary_path('alignments-counts')

        # try to create the temporary directories
        try:
            os.mkdir(temp_alignment_dir)
        except OSError:
            pass
        
        set_bits = self.get_option('set_FLAG_bits')
        unset_bits = self.get_option('unset_FLAG_bits')

        with open(out_alignment_counts_file, 'w') as alignment_counts:
            header = ["SAMPLE","FLAGS_SET", "FLAGS_UNSET", "COUNTS"]
            alignment_counts.write( ",".join(header) + "\n")
        
        for sample, input_path in sample_alignment_paths_dict.iteritems():
            for i in range(len(self.get_option('set_FLAG_bits'))):
                out_file = os.path.join(
                    temp_mapped_dir, "%s-f_%s_-F_%s.txt" % 
                    (sample, set_bits[i], unset_bits[i]))
                
                with process_pool.ProcessPool(self) as pool:
                    with pool.Pipeline(pool) as pipeline:
                        cat4m = [self.get_tool('cat4m'), input_path]
                        samtools = [self.get_tool('samtools'), 'view', '-c',
                                    '-f', unset_bits[i],
                                    '-F', set_bits[i],
                                    '-q', str(self.get_option(
                                        'exclude_MAPQ_smaller_than')),
                                    '-', '-o', '-']
                        
                        pipeline.append(cat4m)
                        pipeline.append(samtools, stdout_path = out_file)
                        
                counts = int()
                with open(out_file) as f:
                    counts = int(f.readline().rstrip())

                with open(out_alignment_counts_file, 'a') as alignment_counts:
                    alignment_counts.write("%s,%s,%s,%s\n" % 
                                           (sample, set_bits[i], 
                                            unset_bits[i], counts))
