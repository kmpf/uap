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
        
        self.require_tool('wc')
        self.require_tool('pigz')
        self.require_tool('cat')

    def declare_runs(self):
        # get a list of all read files we have to count
        sample_input_paths_dict = dict()
        reads_counts_files = dict()
        read_files = list()

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            # Check that all files are fastq.gz files
            for _ in input_paths:
                ext_matches =  self.which_extensions_match_file_path(
                    _, ("fastq", "gz", "gzip"))
                if len(ext_matches) == 0:
                    raise StandardError(
                        "File %s has none of these extensions %s. "
                        "Can't handle it."
                        % (_, ", ".join()) )

            sample_input_paths_dict[run_id] = input_paths
            read_files.extend(input_paths)

        run_id = self.get_step_name()
        with self.declare_run(run_id) as run:
            run.add_output_file('statistics', run_id + 
                                '.raw_read_statistics.txt', read_files)
            run.add_private_info('sample_input_paths', sample_input_paths_dict)
            run.new_exec_group()

    def execute(self, run_id, run):
        out_path = run.get_single_output_file_for_annotation('statistics')
        sample_input_paths_dict = run.get_private_info('sample_input_paths')
        temp_count_dir = self.get_temporary_path('raw-counts')

        try:
            os.mkdir(temp_count_dir)
        except OSError:
            pass

        with open(out_path, 'w') as statistics_file:
            header = ["SAMPLE_NAME", "RAW_READS"]
            statistics_file.write( ",".join(header) + "\n" )
        
            for sample, input_paths in sample_input_paths_dict.iteritems():
                temp_count_file = os.path.join(
                    temp_count_dir, "%s.counts" % sample)

                all_files_gzipped = True
                for _ in input_paths:
                    if len( self.which_extensions_match_file_path(_, ("gz", "gzip")) ) == 0:
                        all_files_gzipped = False
                
                with process_pool.ProcessPool(self) as pool:
                    with pool.Pipeline(pool) as pipeline:
                        if all_files_gzipped:
                            pigz = [self.get_tool('pigz'), '--decompress', 
                                    '--blocksize', '4096', '--processes', 
                                    '2', '-c']
                            pigz.extend(input_paths)
                            pipeline.append(pigz)
                        else:
                            cat = [self.get_tool('cat')]
                            cat.extend(input_paths)
                            pipeline.append(cat)

                        wc = [self.get_tool('wc'), '-l']
                        pipeline.append(wc, stdout_path = temp_count_file)
                
                f = open(temp_count_file)
                counts = (int(f.readline().rstrip()) / 4)
                f.close()

                statistics_file.write( "%s,%s\n" % (sample, counts))
                
