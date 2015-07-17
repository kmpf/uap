import sys
from abstract_step import *
import process_pool
import yaml
import pprint

class Raw_Counts_To_Norm_Counts (AbstractStep):

    def __init__(self, pipeline):
        super(Raw_Counts_To_Norm_Counts, self).__init__(pipeline)

        self.set_cores (10)

        # Input files -- for each sample separately
        self.add_connection('in/counts')
        
        # Output files -- several CSV files combining all samples (raw_counts.csv, raw_counts_median.csv, rfpkms.csv, rfpkms_median.csv, tpms.csv, tpms_median.csv)
        self.add_connection('out/csv_counts_raw')
        self.add_connection('out/csv_counts_raw_median')
        self.add_connection('out/csv_counts_rfpkms')
        self.add_connection('out/csv_counts_rfpkms_median')
        self.add_connection('out/csv_counts_tpms')
        self.add_connection('out/csv_counts_tpms_median')

        # Error and log files
        self.add_connection('out/log_stdout')
        self.add_connection('out/log_stderr')

        # Which tools are required?
        self.require_tool('Rscript')

        # Other options required
        self.add_option('run_id', str)
        self.add_option('path_rscript', str)

    def declare_runs(self):

        # The only information needed from the previous run is the basename of the count files
        basename = []
        count_files = []
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection ('in/counts'):
            basename.append(os.path.basename(input_paths[0]))
            count_files.append(input_paths[0])
            
        run_id = self.get_option('run_id')
        with self.declare_run(run_id) as run:
#            run.add_private_info('basename', basename)
            run.add_output_file('csv_counts_raw', '%s-raw_counts.csv' % run_id, count_files)


    def execute(self, run_id, run):
        
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                
#                 pp = pprint.PrettyPrinter(indent=4)
#                 pp.pprint(basename)

                 ## Call Rscript
                 rscript = [self.get_tool('Rscript'), self.get_option('path_rscript')]
                 
                 #log_stderr = run.get_single_output_file_for_annotation('log_stderr')
                 #log_stdout = run.get_single_output_file_for_annotation('log_stdout')
                 
                 pipeline.append (rscript) #, stderr_path = log_stderr, stdout_path = log_stdout)
