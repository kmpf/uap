import sys
from abstract_step import *
import process_pool
import yaml


class CountMappedReads(AbstractStep):

    def __init__(self, pipeline):
        super(CountMappedReads, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/statistics')
        
        self.require_tool('cat4m')
        self.require_tool('samtools')
        self.require_tool('pigz')

        self.add_option('set_FLAG_bits', list, default=0)
        self.add_option('unset_FLAG_bits', list, default=0)
        self.add_option('exclude_MAPQ_smaller_than', list, default=0)

    def declare_runs(self):
        # set_FLAG_bits and unset_FLAG_bits must have equal length,
        # because both are set everytime samtools is called
        set_bits = self.get_option('set_FLAG_bits')
        unset_bits = self.get_option('unset_FLAG_bits')

        if len(set_bits) != len(unset_bits):
            raise StandardError("set_FLAG_bits has %s elements but " +
                                "unset_FLAG_bits has %s elements. Both " +
                                "options need to have same length.")

        # create runs and check file types
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                # Only one BAM/SAM file expected as input file
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file.")

                basename = os.path.basename(input_paths[0]).split('.')
                if 'sam' in basename:
                    sam_index = basename.index('sam')
                    basename = basename[:sam_index]
                elif 'bam' in basename:
                    bam_index  = basename.index('bam')
                    basename = basename[:bam_index]
                else:
                    raise StandardError("File %s is neither a BAM nor SAM file" % 
                                        (input_paths[0]) )
                basename = '.'.join(basename)
        
                files_dict = {}
                temp_count_files = []
                for i in range(len(set_bits)):
                    count_file = basename + "-f_%s_-F_%s.txt" % (set_bits[i], unset_bits[i])
                    temp_count_files.append(count_file)
                    files_dict[count_file] = {'FLAGS_SET': set_bits[i],
                                              'FLAGS_UNSET': unset_bits[i] }

                run.add_output_file('statistics', basename + 
                                    '.read_statistics.txt', input_paths)
                run.add_private_info('in-alignment', input_paths)                
                run.add_private_info('temp-count-files', temp_count_files)
                run.add_private_info('files-dict', files_dict)

    def execute(self, run_id, run):
        alignment_path = run.get_private_info('in-alignment')
        temp_count_files = run.get_private_info('temp-count-files')
        files_dict = run.get_private_info('files_dict')
        single_count_dir = self.get_temporary_path('single-counts')        

        # samtools view -c 
        # Ich will hier einen Schritt haben in dem ich eine beliebige Anzahl an 
        # Auswertungen machen kann und die dann ordentlich in eine Datei 
        # geschrieben werden!!!

        with process_pool.ProcessPool(self) as pool:
            for i in range(len(set_bits)):
                counts_file = single_count_dir + temp_count_files[i]
                print(counts_file)
                with pool.Pipeline(pool) as pipeline:
                    cat4m = [self.get_tool('cat4m'), alignment_path]
                    samtools = [self.get_tool('samtools'), 'view', '-c',
                                '-f', self.get_option('unset_FLAG_bits')[i],
                                '-F', self.get_option('set_FLAG_bits')[i],
                                '-q', self.get_option('exclude_MAPQ_smaller_than')]
                    samtools.extend(['-', counts_file])
                                        
                    pipeline.append(cat4m)
                    pipeline.append(samtools)

        # Read in all count files and create the final statistics file
        header = ["FLAGS_SET", "FLAGS_UNSET"]
        statistics_file = open(
            run.get_single_output_file_for_annotation('statistics'), 'w')
        statistics_file.write( header.join(",") + ",COUNTS")
        for count_file in files_dict:
            f = open(single_count_dir + count_file)
            counts = f.readline()
            line = "%s, %s, %s"  % ( files_dict[count_file][header[0]],
                                     files_dict[count_file][header[1]], counts)
            statistics_file.write( line )
