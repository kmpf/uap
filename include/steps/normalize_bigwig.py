import sys
import yaml

from ..abstract_step import *
from .. import process_pool

### TODO
### Each run_id is associated with two wig files plus and minus
### Need to take care of both

class Normalize_BigWig (AbstractStep):
    '''
    Normalize_wig from the RSeQC package normalizes wig files based on wigsum.

    From the RSeqQC package documentation (ver 3.2.9):
    Therefore, only normalized to "total read count" is problematic if read length is different between samples.
    Here we normalize every bigwig file into the same wigsum. wigsum is the summary of signal value across the genome. for example, wigsum = 100,000,000 equals to
    the coverage achieved by 1 million 100nt long reads or 2 million 50nt long reads.

    Options:
    --version 	show program version number and exit
    -h, --help 	show this help message and exit
    -i BIGWIG_FILE, --bwfile=BIGWIG_FILE
     	Input BigWig file. [required]
    -o OUTPUT_WIG, --output=OUTPUT_WIG
     	Output wig file. [required]
    -s CHROMSIZE, --chromSize=CHROMSIZE
     	Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name, second column is size of the chromosome. [required]
    -t TOTAL_WIGSUM, --wigsum=TOTAL_WIGSUM
     	Specified wigsum. 100000000 equals to coverage of 1 million 100nt reads. default=100000000 [optional]
    -r REFGENE_BED, --refgene=REFGENE_BED
     	Reference gene model in bed format. [optional]
    -c CHUNK_SIZE, --chunk=CHUNK_SIZE
     	Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=100000 (bp) [optional]

    All options are available through this interface.
    The chromsizes file can be obtain using UCSC tool fetchChromSizes.
    -i and -o do not apply as they are set by the pipeline automatically
    '''

    def __init__(self, pipeline):
        super(Normalize_BigWig, self).__init__(pipeline)

        self.set_cores (1)

        self.add_connection('in/tracks')
        self.add_connection('out/tracks')
        

        self.add_option('s', str, optional=False)
        self.add_option('t', int, optional=True)
        self.add_option('r', str, optional=True)
        self.add_option('c', int, optional=True)

        self.require_tool('normalize_bigwig')
        self.require_tool('wigToBigWig')
        

    def declare_runs (self):
      ## check availability of chromsizes files
      chromSizes = self.get_option ('s')
      if not os.path.exists(chromSizes):
        raise StandardError ('Cannot find chromSizesFile %s' % chromSizes)
    
      ## if r is set, check file availability of r
      if self.is_option_set_in_config('r'):
          refBed=self.get_option('r')
          if not os.path.exists(refBed):
              raise StandardError ('Cannot find reference bed file  %s' % refBed)
          
      for run_id, input_paths in self.get_run_ids_and_input_files_for_connection ('in/tracks'):
          ## Check whether we have one ot more bigwig files for a run_id,
          ## usually we have two and need to declare two run_ids
          if (len(input_paths) > 1):
              
              ## Declare two runs one for the plus one for the minus bigwig file
              ## use the path without .bw as runid
              for path in input_paths:
                  new_run_id=os.path.basename(path)
                  new_run_id=new_run_id.replace(".bw", "")
                  new_run_id=new_run_id.replace(".","-")
                  
                  with self.declare_run (new_run_id) as run:
                      run.add_output_file('tracks', '%s-norm.bw'  % new_run_id, [path])
                      run.new_exec_group()
          else:
              with self.declare_run(run_id) as run:
                  run.add_output_file('tracks', '%s-norm.bw'  %
                                      run_id, input_paths)
                  run.new_exec_group()

          #print new_run_id
          #print run.get_single_output_file_for_annotation ('tracks')
          #norm_bw_path = run.get_single_output_file_for_annotation ('tracks')
          #print  run.get_input_files_for_output_file(norm_bw_path)
          

    def execute (self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                ## Obtain the output files
                ## RSeQC normalize_wig produces a wig file not bigwig as output need to generate bigwig subsequently
                norm_wig_path=self.get_temporary_path('norm_wig', 'output')
                norm_bw_path = run.get_single_output_file_for_annotation ('tracks')

                ## Obtain input file for output
                in_bw = run.get_input_files_for_output_file(norm_bw_path)[0]

                #print ("in_bw:", in_bw)
                
                
                #sorted_bam_path = run.get_single_output_file_for_annotation('alignments')
                #sorted_bai_path = run.get_single_output_file_for_annotation('indices')
                #unsorted_bam_path = self.get_temporary_path('merge_bam_unsorted', 'output')
                
                normalize_wig = [self.get_tool('normalize_bigwig')]
                normalize_wig.extend (['-i', in_bw,
                                       '-o', norm_wig_path,
                                       '-s', self.get_option('s')])
                
                if self.is_option_set_in_config('t'):
                    normalize_wig.extend  (['-t', self.get_option('t')])

                if self.is_option_set_in_config('r'):
                    normalize_wig.extend  (['-r', self.get_option('r')])

                if self.is_option_set_in_config('c'):
                    normalize_wig.extend  (['-c', self.get_option('c')])

                #def warning(*objs):
                #    print("WARNING: ", *objs, file=sys.stderr)
                #warning(normalize_wig)
                
                pipeline.append (normalize_wig)
                

               
        ## Convert normalized wig file to bigwig
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                wigToBigWig = [self.get_tool('wigToBigWig'), norm_wig_path, self.get_option('s'), norm_bw_path]
                
                
                pipeline.append(wigToBigWig)


        if os.path.exists(norm_wig_path):
            os.unlink(norm_wig_path)
