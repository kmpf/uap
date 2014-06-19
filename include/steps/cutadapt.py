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
        
        self.add_connection('in/reads')
        self.add_connection('out/reads')
        self.add_connection('out/log')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('cutadapt')
        self.require_tool('fix_qnames')
        
        self.add_option('adapter-type', str, optional=True)
        self.add_option('adapter-R1', str, optional = False)
        self.add_option('adapter-R2', str, optional = True)
#        self.add_option('adapter', str, optional = True)
        self.add_option('fix_qnames', bool, default = False)
        
    def declare_runs(self):
        # fetch all incoming run IDs which produce reads...
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            is_paired_end = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')

            ## Make sure the adapter type is one of a, b or g 
            if self.is_option_set_in_config('adapter-type'):
                if not (self.get_option('adapter-type') in set(['a','b','g'])):
                    raise StandardError("Option 'adapter-type' must be either 'a','b', or 'g'!")
                                
            # make sure that adapter-R1/adapter-R2 are set correctly
            # according to paired_end info... this kind of mutual exclusive option
            # checking is quite complicated, so we do it here.

            if is_paired_end:
                if not self.is_option_set_in_config('adapter-R2'):
                    raise StandardError("Option 'adapter-R2' required because sample %s is paired end!" % run_id)
            elif self.is_option_set_in_config('adapter-R2'):
                    raise StandardError("Option 'adapter-R2' not allowed because sample %s is not paired end!" % run_id)

            # decide which read type we'll handle based on whether this is
            # paired end or not
            read_types = list()
            if not is_paired_end:
                read_types = ['-R1']
            else:
                read_types = ['-R1', '-R2']

            # put input files into R1/R2 bins (or one single R1 bin)
            input_path_bins = dict()
            for _ in read_types:
                input_path_bins[_] = list()
            for path in input_paths:
                which = '-R1'
                if is_paired_end:
                    which = '-' + misc.assign_string(os.path.basename(path), ['R1', 'R2'])
                input_path_bins[which].append(path)

            # now declare two runs
            for which in read_types:
                with self.declare_run("%s%s" % (run_id, which)) as run:
                    # add paired end information
                    run.add_private_info('paired_end', is_paired_end)
                    if is_paired_end:
                        run.add_private_info('paired_end_read', which.replace('-', ''))

                    # add adapter information, insert correct index first if necessary
                    adapter = self.get_option('adapter%s' % which)
                    
                    if '((INDEX))' in adapter:
                        index = self.find_upstream_info_for_input_paths(input_paths, 'index%s' % which)
                        if which == '-R2':
                            complements = string.maketrans('acgtACGT', 'tgcaTGCA')
                            index = index.translate(complements)[::-1]

                        adapter = adapter.replace('((INDEX))', index)
                     # make sure the adapter is looking good
                    if re.search('^[ACGT]+$', adapter) == None:
                        raise StandardError("Unable to come up with a legit-looking adapter: " + adapter)
                    run.add_private_info('adapter', adapter)

                    # add output files: reads and log
                    run.add_output_file("reads", "%s-cutadapt%s.fastq.gz" % (run_id, which), input_path_bins[which])
                    run.add_output_file("log", "%s-cutadapt%s-log.txt" % (run_id, which), input_path_bins[which])


    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                out_path = run.get_single_output_file_for_annotation('reads')

                # set up processes
                cat4m = [self.get_tool('cat4m')]
                cat4m.extend(sorted(run.get_input_files_for_output_file(out_path)))
                if self.is_option_set_in_config('adapter-type'):
                    adapterType = '-' + self.get_option('adapter-type')
                else:
                    adapterType='-a'
                
                pigz1 = [self.get_tool('pigz'), '--processes', '1', '--decompress', '--stdout']
                
                fix_qnames = [self.get_tool('fix_qnames')]

                cutadapt = [self.get_tool('cutadapt'), adapterType, run.get_private_info('adapter'), '-']

                pigz2 = [self.get_tool('pigz'), '--blocksize', '4096', '--processes', '1', '--stdout']

                # create the pipeline and run it
                pipeline.append(cat4m)
                pipeline.append(pigz1)
                if self.get_option('fix_qnames') == True:
                    pipeline.append(fix_qnames)
                pipeline.append(cutadapt, stderr_path = run.get_single_output_file_for_annotation('log'))
                pipeline.append(pigz2, stdout_path = out_path)
