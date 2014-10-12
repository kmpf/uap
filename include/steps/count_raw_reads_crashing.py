import sys
from abstract_step import *
import process_pool
import yaml


class CountRawReads(AbstractStep):

    def __init__(self, pipeline):
        super(CountRawReads, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/reads')
        self.add_connection('out/statistics')
        
        self.require_tool('wc')
        self.require_tool('pigz')

    def declare_runs(self):
        # get a list of all read files we have to count
        sample_files_dict = dict()
        reads_counts_files = dict()
        read_files = list()
        # Check that all files are fastq.gz files
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            sample_files_dict[run_id] = input_paths
            read_files.extend(input_paths)
            for f in input_paths:
                if not f.endswith("fastq.gz"):
                    raise StandardError("Input file %s is not ending with "
                                        "'fastq.gz'." % f)
                basename = os.path.basename(f)
                reads_counts_files[f] = basename.replace("fastq.gz", "counts")

        # Create new run ID name
        run_id = "%s" % ( self.get_step_name() )
        print(run_id)

        # Create the new single run (just one)
        with self.declare_run(run_id) as run:
            # Define the output file which will be created
            print(len(read_files))
            run.add_output_file('statistics', run_id + 
                                '.raw_read_statistics.txt', read_files)
            # Store the map of run_id to input_paths as private info
            run.add_private_info('read-files', sample_files_dict)
            # Store the map of input_files to count_files as private info
            run.add_private_info('reads-counts-files', reads_counts_files)

    def execute(self, run_id, run):
        print("Over here!")
        read_files = run.get_private_info('read-files')
        reads_counts_files = run.get_private_info('reads-counts-files')

        # Initialize process pool
        with process_pool.ProcessPool(self) as pool:
            # Each file from each sample is counted and the result 
            # (a single number) is stored in a temporary file
            for sample, input_paths in read_files:
                for raw_reads in input_paths:
                    with pool.Pipeline(pool) as pipeline:
                
                        pigz = [self.get_tool('pigz'), '--decompress', 
                                '--blocksize', '4096', '--processes', '2', '-c']
                        pipeline.append(pigz)
                        echo = [self.get_tool('echo'), '$((`wc -l`/4))']
                        pipeline.append(
                            echo, stdout_path = reads_counts_files[raw_reads])

        # Read in all count files and create the final statistics file
        header = ["SAMPLE_NAME", "RAW_READS"]
        statistics_file = open(
            run.get_single_output_file_for_annotation('statistics'), 'w')
        statistics_file.write( header.join(",") )
        # Collect all results for each sample, add them up and store them in the
        # final result file
        results = dict()
        for sample, input_paths in read_files:
            results = 0
            for raw_reads in input_paths:
                f = open(reads_counts_files[raw_reads])
                counts = f.readline()
                results += counts
            statistics_file.write( "%s, %s" % (sample, results) )
