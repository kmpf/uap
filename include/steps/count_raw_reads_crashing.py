import sys
from abstract_step import *
import process_pool
import yaml

# Taken from:
# http://stackoverflow.com/questions/2023608/check-what-files-are-open-in-python
#import __builtin__
#openfiles = set()
#oldfile = __builtin__.file
#class newfile(oldfile):
#    def __init__(self, *args):
#        self.x = args[0]
#        print "### OPENING %s ###" % str(self.x)
#        oldfile.__init__(self, *args)
#        openfiles.add(self)
#        printOpenFiles()
#
#    def close(self):
#        print "### CLOSING %s ###" % str(self.x)
#        oldfile.close(self)
#        openfiles.remove(self)
#oldopen = __builtin__.open
#def newopen(*args):
#    return newfile(*args)
#__builtin__.file = newfile
#__builtin__.open = newopen
#
#def printOpenFiles():
#    print "### %d OPEN FILES: [%s]" % (len(openfiles), ", ".join(f.x for f in openfiles))
#    print "### %d OPEN FILES" % len(openfiles)

class CountRawReads(AbstractStep):

    def __init__(self, pipeline):
        super(CountRawReads, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/reads')
        self.add_connection('out/statistics')
        
        self.require_tool('echo')
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
 
        # Create the new single run (just one)
        with self.declare_run(run_id) as run:
            # Define the output file which will be created
            run.add_output_file('statistics', run_id + 
                                '.raw_read_statistics.txt', read_files)
            # Store the map of run_id to input_paths as private info
            run.add_private_info('read-files', sample_files_dict)
            # Store the map of input_files to count_files as private info
            run.add_private_info('reads-counts-files', reads_counts_files)

    def execute(self, run_id, run):
        read_files = run.get_private_info('read-files')
        reads_counts_files = run.get_private_info('reads-counts-files')
        temp_count_dir = self.get_temporary_path('raw-counts')

        try:
            os.mkdir(temp_count_dir)
        except OSError:
            pass

        statistics_file = open(
            run.get_single_output_file_for_annotation('statistics'), 'w')
        header = ["SAMPLE_NAME", "RAW_READS"]
        statistics_file.write( ",".join(header) + "\n" )
        statistics_file.close()
        
        counts_per_sample = int
        for sample, input_paths in read_files.iteritems():
            # Reset counter to zero for every sample
            counts_per_sample = 0
            for raw_reads in input_paths:
                # Initialize process pool
                with process_pool.ProcessPool(self) as pool:
                    temp_count_file = os.path.join(
                        temp_count_dir, reads_counts_files[raw_reads])
                    print(temp_count_file)
                    # Each file from each sample is counted and the result 
                    # (a single number) is stored in a temporary file
                    with pool.Pipeline(pool) as pipeline:
                        
                        pigz = [self.get_tool('pigz'), '--decompress', 
                                '--blocksize', '4096', '--processes', '2', '-c',
                                raw_reads]
                        wc = [self.get_tool('wc'), '-l']
                        
                        pipeline.append(pigz)
                        pipeline.append(wc, stdout_path = temp_count_file)

                    # Collect all results for each sample, add them up and store 
                    # them in the final result file

                f = open(temp_count_file)
                counts = int(f.readline().rstrip())
                f.close()
                counts_per_sample += (counts / 4)

            print(sample)
            statistics_file = open(
                run.get_single_output_file_for_annotation('statistics'), 'w')
            statistics_file.write( "%s, %s\n" % (sample, counts_per_sample) )
            statistics_file.close()
