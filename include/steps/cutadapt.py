import sys
from abstract_step import *
import pipeline
import re
import process_pool
import yaml

class Cutadapt(AbstractStep):
    '''
    The cutadapt step can be used to clip adapter sequences from RNASeq reads.
    
    Any adapter may contain ``((INDEX))`` which will be replaced with every
    sample's index. The resulting adapter is checked for sanity and an
    exception is thrown if the adapter looks non-legit.
    '''
    
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)
        
        self.set_cores(3)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        self.add_connection('out/log_first_read')
        self.add_connection('out/log_second_read')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('cutadapt')
        self.require_tool('fix_qnames')

        # Options for cutadapt
        self.add_option('adapter-type', str, optional=True)
        self.add_option('adapter-R1', str, optional = True)
        self.add_option('adapter-R2', str, optional = True)
        self.add_option('adapter-file', str, optional=True)
        self.add_option('use_reverse_complement', bool, default=False)
        self.add_option('minimal-length', int, default=10, optional=True)

        self.add_option('fix_qnames', bool, default = False)
        
    def declare_runs(self):
        found_files = dict()
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        paired_end_info = dict()
        # fetch all incoming run IDs which produce reads...
        for read in read_types.keys():
            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/%s' % read):
                if input_paths != [None]:
                    paired_end_info[run_id] = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
                    # make sure that adapter-R1/adapter-R2 or adapter-file are set
                    # correctly according to paired_end info... this kind of mutual 
                    # exclusive option checking is a bit tedious, so we do it here.
                    
                    if paired_end_info[run_id]:
                        if ( not self.is_option_set_in_config('adapter-R2') and 
                             not self.is_option_set_in_config('adapter-file') ):
                            raise StandardError(
                                "Option 'adapter-R2' or " +
                                "'adapter-file' required because " +
                                "sample %s is paired end!" % run_id)

                    elif ( self.is_option_set_in_config('adapter-R2') and
                           not self.is_option_set_in_config('adapter-file') ):
                        raise StandardError(
                            "Option 'adapter-R2' not allowed because " +
                            "sample %s is not paired end!" % run_id)
                        
                    if ( self.is_option_set_in_config('adapter-file') and
                         self.is_option_set_in_config('adapter-R1') ):
                        raise StandardError(
                            "Option 'adapter-R1' and 'adapter-file' " +
                            "are both set but are mutually exclusive!" )

                    if ( not self.is_option_set_in_config('adapter-file') and
                         not self.is_option_set_in_config('adapter-R1') ):
                        raise StandardError(
                            "Option 'adapter-R1' or 'adapter-file' " +
                            "required to call cutadapt for sample %s!" % 
                            run_id)



                    # save information per file in found_files
                    if not run_id in found_files:
                        found_files[run_id] = dict()
    
                    if not read in found_files[run_id]:
                        found_files[run_id][read] = list()
                    found_files[run_id][read].extend(input_paths)
    
                    ## Make sure the adapter type is one of a, b or g 
                    if self.is_option_set_in_config('adapter-type'):
                        if not (self.get_option('adapter-type') in set(['a','b','g'])):
                            raise StandardError("Option 'adapter-type' must be "
                                "either 'a','b', or 'g'!")
                

        # now declare two runs
        for run_id in found_files.keys():
            with self.declare_run(run_id) as run:
                run.new_exec_group()
                # add paired end information
                run.add_private_info('paired_end', paired_end_info[run_id])
                for read in found_files[run_id].keys():
                    # find correct adapter information and add info to run
                    if self.is_option_set_in_config('adapter%s' % read_types[read]):
                        # add adapter information, insert correct index first if 
                        # necessary
                        adapter = self.get_option('adapter%s' % read_types[read])
                        
                        # add index to adapter sequence if necessary
                        if '((INDEX))' in adapter:
                            index = self.find_upstream_info_for_input_paths(
                                found_files[run_id][read], 
                                'index%s' % read_types[read])
                            adapter = adapter.replace('((INDEX))', index)
                        
                        # create reverse complement if we are asked for it
                        if self.get_option('use_reverse_complement'):
                            complements = string.maketrans('acgtACGT', 'tgcaTGCA')
                            adapter = adapter.translate(complements)[::-1]
                            
                        # make sure the adapter is looking good
                        if re.search('^[ACGT]+$', adapter) == None:
                            raise StandardError("Unable to come up with a "+
                                                "legit-looking adapter: " + 
                                                adapter)
                        run.add_private_info('adapter-%s' % read, adapter)
                    elif self.is_option_set_in_config('adapter-file'):
                        adapter = "file:" + self.get_option('adapter-file')
                        run.add_private_info('adapter-%s' % read, adapter)

                    # add output files: reads and log
                    run.add_output_file(read, 
                        "%s-cutadapt%s.fastq.gz" % (run_id, read_types[read]), 
                        found_files[run_id][read])
                    run.add_output_file("log_%s" % read, 
                        "%s-cutadapt%s-log.txt" % (run_id, read_types[read]), 
                        found_files[run_id][read])

    def execute(self, run_id, run):
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        for read in read_types.keys():
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    out_path = run.get_single_output_file_for_annotation(read)

                    # set up processes
                    # cat all files
                    cat4m = [self.get_tool('cat4m')]
                    cat4m.extend(sorted(run.get_input_files_for_output_file(out_path)))
                    pipeline.append(cat4m)
                        
                    # Decompress the input files
                    pigz1 = [self.get_tool('pigz'), '--processes', '1', 
                             '--decompress', '--stdout']
                    pipeline.append(pigz1)
                    
                    # fix qnames if wanted
                    fix_qnames = [self.get_tool('fix_qnames')]
                    if self.get_option('fix_qnames') == True:
                        pipeline.append(fix_qnames)

                    # set option for adapter clipping
                    if self.is_option_set_in_config('adapter-type'):
                        adapterType = '-' + self.get_option('adapter-type')
                    else:
                        adapterType='-a'

                    # Clip adapters from piped data
                    cutadapt = [self.get_tool('cutadapt'), adapterType, 
                                run.get_private_info('adapter-%s' % read), '-']
                    pipeline.append(cutadapt, stderr_path = 
                                    run.get_single_output_file_for_annotation(
                                        'log_%s' % read))
                    
                    # Compress output to file
                    pigz2 = [self.get_tool('pigz'), '--blocksize', '4096', 
                             '--processes', '1', '--stdout']
                    pipeline.append(pigz2, stdout_path = out_path)
