import sys
from abstract_step import *
import process_pool
import yaml
import os
from logging import getLogger

logger=getLogger('uap_logger')

class sra_fastq_dump (AbstractStep):
    '''
    sra tools is a suite from NCBI to handle sra (short read archive) files.
    fastq-dump is an sra tool that dumps the content of an sra file in fastq format

    The following options cannot be set, as they would interefere with the pipeline implemented in this step
     -O|--outdir <path>               Output directory, default is working 
                                   directory '.' ) 
     -Z|--stdout                      Output to stdout, all split data become 
                                   joined into single stream 
     --gzip                           Compress output using gzip 
     --bzip2                          Compress output using bzip2 

      Multiple File Options              Setting these options will produce more
                                     than 1 file, each of which will be suffixed
                                     according to splitting criteria.
       --split-files                    Dump each read into separate file.Files 
                                   will receive suffix corresponding to read 
                                   number 
       --split-3                        Legacy 3-file splitting for mate-pairs: 
                                        First biological reads satisfying dumping 
                                        conditions are placed in files *_1.fastq and 
                                        *_2.fastq If only one biological read is 
                                        present it is placed in *.fastq Biological 
                                        reads and above are ignored. 
       -G|--spot-group                  Split into files by SPOT_GROUP (member name) 
       -R|--read-filter <[filter]>      Split into files by READ_FILTER value 
                                        optionally filter by value: 
                                        pass|reject|criteria|redacted 
       -T|--group-in-dirs               Split into subdirectories instead of files 
       -K|--keep-empty-files            Do not delete empty files 


    Details on fastq-dump can be found at
    https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

    To make IO cluster friendly, fastq-dump is not reading th sra file directly. Rather, dd with configurable blocksize is used
    to provide the sra file via a fifo to fastq-dump- 

    The executed calls lools like this 
    
    mkfifo sra_fifo
    dd bs=4M if=<sra-file> of=sra_fifo
    fastq-dump -Z sra_fifo | pigz --blocksize 4096 --processes 2 > file.fastq

    '''

    def __init__(self, pipeline):
    super (sra_fastq_dump, self).__init__(pipeline)
    self.set_cores(10) # set # of cores for cluster, it is ignored if run locally

    # connections - indentifier for in/output
    self.add_connection('in/sequence')
    self.add_connection('out/sequence')

    self.require_tool('fastq-dump')
    self.require_tool('dd')
    self.require_tool('pigz')
    self.require_tool('mkfifo')

    

   ## Complete set of options (excluding those mentioned above) for version 2.7.0
   ## INPUT
    self.add_option('accession', str, optional=True,
                    description="Replaces accession derived from <path> in filename(s) and deflines (only for single table dump)")
    
    self.add_option('table', str, optional=True,
                    description="Table name within cSRA object, default is \"SEQUENCE\"")

    ## Processing
    ### Read Splitting
    self.add_option('split-spot', str, optional=True,
                    description="Split spots into individual reads")
    ### Full Spot Filters
    self.add_option('minSpotId', int, optional=True,
                    description="Minimum spot id")

    self.add_option('maxSpotId', int, optional=True,
                    description="Maximum spot id")

    self.add_option('spot-groups', str, optional=True,
                    description="Filter by SPOT_GROUP (member): name[,...]")

    self.add_option('clip', bool, optional=True,
                    description='Apply left and right clips')

    ### Common Filters
    self.add_option('minReadLen', int, optional=True,
                    description='Filter by sequence length >= <len>')

    self.add_option('read-filter', str, optional=True,
                    description='Split into files by READ_FILTER value optionally filter by value: pass|reject|criteria|redacted')

    self.add_option('qual-filter', bool, optional=True,
                    description='Filter used in early 1000 Genomes data: no sequences starting or ending with >= 10N')

    self.add_option('qual-filter-1', bool, optional=True,
                    description='Filter used in current 1000 Genomes data')

    ### Filters based on alignments
    self.add_option('aligned', bool, optional=True,
                    description='Dump only aligned sequences')

    self.add_option('unaligned', bool, optional=True,
                    description='Dump only unaligned sequences')

    self.add_option('aligned-region', str, optional=True,
                    description='Filter by position on genome. Name can either be accession.version (ex:NC_000001.10) or file specific name (ex:"chr1" or "1"). "from" and "to" are 1-based coordinates. <name[:from-to]>')

    self.add_option('matepair-distance' str, optional=True,
                    description='Filter by distance beiween matepairs. Use "unknown" to find matepairs split between the references. Use from-to to limit matepair distance on the same reference. <from-to|unknown>')

    ### Filters for individual reads
    self.add_option('skip-technical', bool, optional=True,
                    description='Dump only biological reads')

    ## FORMATTING
    ### Sequence
    self.add_option('dumpcs', bool, optional=True,
                    description='Formats sequence using color space (default for SOLiD),\"cskey\" may be specified for translation.'
    )

    self.add_option('dumpbase', bool, optional=True,
                    description='Formats sequence using base space (default for other than SOLiD).')

    ### Quality
    self.add_option('offset', int, optional=True,
                    description='Offset to use for quality conversion, default is 33')

    self.add_option('fasta', int, optional=True,
                    description='FASTA only, no qualities, optional line wrap width (set to zero for no wrapping). <[line width]>')

    self.add_option('suppress-qual-for-cskey', bool, optional=True,
                    description='supress quality-value for cskey ')

    ### Defline
    self.add_option('origfmt', bool, optional=True,
                    description='Defline contains only original sequence name')

    self.add_option('readids', bool, optional=True,
                    description='Append read id after spot id as \'accession.spot.readid\' on defline.')

    self.add_option('helicos', bool, optional=True,
                    description='Helicos style defline')

    self.add_option('defline-seq', str, optional=True,
                    description='Defline format specification for sequence.  <fmt> is string of characters and/or variables. The variables can be one of: $ac - accession, $si spot id, $sn spot name, $sg spot group (barcode), $sl spot length in bases, $ri read number, $rn read name, $rl read length in bases. '[]' could be used for an optional output: if all vars in [] yield empty values whole group is not printed. Empty value is empty string or for numeric variables. Ex: @$sn[_$rn]/$ri \'_$rn\' is omitted if name is empty. <fmt>')

    self.add_option('defline-qual', str, optional=True,
                    description='Defline format specification for quality.')

    
    ## OTHER
    self.add_option('disable-multithreading', bool, optional=True,
                    description='disable multithreading')

    self.add_option('log-level', str,  optional=True,
                    description='Logging level as number or enum string One of (fatal|sys|int|err|warn|info) or (0-5). Current/default is warn. <level>')

    self.add_option('verbose', bool, optional=True,
                    description='Increase the verbosity level of the program. Use multiple times for more verbosity'.)

    self.add_option('ncbi_error_report', str, optional=True,
                    description='Control program execution environment report generation (if implemented). One of (never|error|always). Default is error. <error>')

    self.add_option('legacy-report', bool, optional=True,
                    description='use legacy style 'Written spots' for tool')

    # [Options for 'dd':]
    self.add_option('dd-blocksize', str, optional = True, default = "256k")
    self.add_option('max_cores', int, optional = True, default= 10, description='Maximum number of cores available on the cluster')
    #########

    def runs (self, run_ids_connections_files):

        options = ['accession','table','split-spot','minSpotId','maxSpotId', 'spot-groups',
                   'clip','minReadLen', 'read-filter','qual-filter','aligned','unaligned',
                   'aligned-region', 'matepair-distance', 'skip-technical', 'dumpcs',
                   'dumpbase','offset','fasta','suppress-qual-for-cskey', 'origfmt',
                   'readids', 'helicos', 'defline-seq', 'defline-qual','disable-multithreading',
                   'log-level','verbose','ncbi_error_report', 'legacy-report', 'dd-blocksize',
                   'max_cores']

        set_options = [option for option in options if \
               self.is_option_set_in_config(option)]
        
        option_list = list()

        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
                else:
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:

                input_paths = run_ids_connections_files[run_id]['in/sequence']
                # check, if only a single input file is provided
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one sra file., but got this %s" % input_paths)

                with run.new_exec_group() as exec_group
                    # 1. Create FIFO for reading sra file
                    fifo_path_sra=run.add_temporary_file (
                        'sra-fifo', designation='input')
                    mkfifo_sra= [self.get_tool('mkfifo'), fifo_path_sra]
                    exec_group.add_command(mkfifo_sra)

                    # 2. Read sra file and output to FIFO
                    dd_sra=[self.get_tool('dd'),
                            'bs=%s' % self.get_option('dd-blocksize'),
                            'if=%s' % input_paths[0],
                            'of=%s' % fifo_path_sra]
                    exec_group.add_command(dd_sra)

                    # 3. Run fastq-dump

                    with exec_group.add_pipeline() as fastq_dump_pipe:
                    fastq_dump=[self.get_tool('fastq-dump'), '--stdout']
                    fastq_dump.extend(option_list)
                    fastq_dump.extend(fifo_path_sra)
                    
                    fastq_dump_pipe.add_command (
                        fastq_dump,
                        stderr_path = run.add_output_file(
                            'log',
                            '%s-fastq-dump-log.txt' % run_id,
                            input_paths)
                    )

                    # 4. Run pigz to efficiently compress fastq file
                    pigz_sequence = [
                        self.get_tool('pigz'),
                        '--stdout',
                        '--blocksize%s' % self.get_option('dd-blocksize'),
                        '--processes', 2 if (self.get_option('max_cores') >=2) else self.get_option('max_cores')
                        ]
                    
                    fastq_dump_pipe.add_command (
                        pigz_sequence,
                        stdout_path = run.add_output_file(
                                'sequence',
                                '%s.fastq.gz' % run_id,
                                input_paths)
                        )

                
                    

                    
