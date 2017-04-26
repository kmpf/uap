import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Segemehl2017(AbstractStep):
    '''
    segemehl_2017 is a software to map short sequencer reads to reference genomes. 
    Unlike other methods, segemehl is able to detect not only mismatches but 
    also insertions and deletions. Furthermore, segemehl is not limited to a 
    specific read length and is able to mapprimer- or polyadenylation 
    contaminated reads correctly.

    This step is a wrapper for an unpublished version of segemehl that we re-
    ceived from Steve Hoffmann on April, 24th 2017. 
    It automatically generates the bed file containing the realigned split 
    reads (split junctions). 
    No need to run testrealign.x after segemehl anymore.

    This step creates at first two FIFOs. The first is used to provide the
    genome data for segemehl and the second is used for the output of the
    unmapped reads::

       mkfifo genome_fifo unmapped_fifo
       cat <genome-fasta> -o genome_fifo

    The executed segemehl command is this::

        segemehl -d genome_fifo -i <genome-index-file> -q <read1-fastq> 
                 [-p <read2-fastq>] -u unmapped_fifo -H 1 -t 11 -s -S -D 0
                 -o /dev/stdout |  pigz --blocksize 4096 --processes 2 -c

    The unmapped reads are saved via these commands::

        cat unmapped_fifo | pigz --blocksize 4096 --processes 2 -c > 
        <unmapped-fastq>

    '''

    def __init__(self, pipeline):
        super(Segemehl2017, self).__init__(pipeline)

        self.set_cores(10) # set # of cores for cluster, it is ignored if run locally

        # connections - indentifier for in/output
        #             - expects list, maybe empty or 'none', e.g. if only first_read info available
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/log')

        # tools from tools section of YAML config file
        # -> logs version information
        # -> for module load/unload stuff
        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('fix_qnames')
        self.require_tool('mkfifo')
        self.require_tool('pigz')
        self.require_tool('segemehl')

        # Options for additional programs
        # these options can be used in YAML config file
        self.add_option('fix-qnames', bool, optional=True, default=False,
                        description="The QNAMES field of the input will "
                        "be purged from spaces and everything thereafter.")

        # Options to set segemehl flags
        ## [INPUT]
        self.add_option('database', str, optional=False,
                        description = "(Space separated list of ) filename(s) of database (e.g. genome) sequene(s)")
#        self.add_option('genome', str, optional=False,
#                        description="Path to genome file")
        # query == input fastq files for R1
        # mate == input fastq files for R2 (if paired end sequencing!)

        self.add_option('index', str, optional=True, default = "none",
                        description="Path to database index for segemehl (default:none)")
        self.add_option('index2', str, optional=True, default = "none",
                        description="Path to second database index for segemehl (default:none)")
        self.add_option('filebins', str, optional=True, default = "none",
                        description = "file bins with basename <string> for easier "
                        "data handling (default:none)")
        self.add_option('bisulfite', int, choices=[0, 1, 2], optional=True, default = 0,
                        description="bisulfite mapping with methylC-seq/Lister "
                        "et al. (=1) or bs-seq/Cokus et al. protocol (=2) "
                        "(default:0)")
        ## [GENERAL]
        self.add_option('minsize', int, optional=True, default = 12,
                        description="minimum size of queries (default:12)")
        # progressbar - we ommit this, since it will only generate larger log files
        self.add_option('brief', bool, optional = True, default = False,
                        description = "brief output")
        self.add_option('checkidx', bool, optional = True, default = False,
                        description = "check index")
        self.add_option('briefcigar', bool, optional = True, default = False,
                        description = "brief cigar string (M vs X and =)")
        self.add_option('threads', int, default=1, optional=True,
                        description="start <n> threads (default:1)")
        # outfile => take care of that later!
        # nomatchfilename
        self.add_option('readgroupfile', str, optional=True, default = "none",
                        description = "filename to read @RG header (default:none)")
        self.add_option('readgroupid', str, optional=True, default = "none",
                        description = "read group id (default:none)")
        
        ## [SEEDPARAMS]
        self.add_option('differences', int, default=1, optional=True,
                        description="search seeds initially with <n> "
                        "differences (default:1)")
        self.add_option('jump', int, optional=True, default = 0, 
                        description= "search seeds with jump size <n> (0=automatic) "
                        "(default:0)")
        self.add_option('nosuflinks', bool, optional = True,
                        description = "dont use suflinks (does not affect index construction, "
                        "for short reads only, increases runtime!)")
        self.add_option('evalue', float, optional=True, 
                        description= "max evalue (default:5.000000)")
        self.add_option('maxsplitevalue', float, optional=True, 
                        description = "max evalue for splits (default:50.000000)")
        self.add_option('maxinterval', int, optional=True, description=
                        "maximum width of a suffix array interval, i.e. a query "
                        "seed will be omitted if it matches more than <n> times "
                        "(default:100)")
        self.add_option('maxout', int, optional=True,
                        description = "maximum number of alignments that will be reported. "
                        "If set to zero, all alignments will be reported (default:0)")
        self.add_option('splits', bool, default=None, optional=True,
                        description="detect split/spliced reads (default:none)")
        self.add_option('MEOP', bool, optional=True, description=
                        "output MEOP field for easier variance calling in SAM "
                        "(XE:Z:)")
        self.add_option('nohead', bool, optional=True, 
                        description = "do not output header")

        ## [SEEDEXTENSIONPARAMS]
        self.add_option('extensionpenalty', int, optional=True, description=
                        "penalty for a mismatch during extension (default:4)")
        self.add_option('dropoff', int, optional=True, description=
                        "dropoff parameter for extension (default:8)")

        ## [ALIGNPARAMS]
        self.add_option('accuracy', int, optional=True, description=
                        "min percentage of matches per read in semi-global "
                        "alignment (default:90)")
        self.add_option('minsplicecover', int, optional=True, description=
                        "min coverage for spliced transcripts (default:80)")
        self.add_option('minfragscore', int, optional=True, description=
                        "min score of a spliced fragment (default:18)")
        self.add_option('minfraglen', int, optional=True, description=
                        "min length of a spliced fragment (default:20)")
        self.add_option('splicescorescale', float, optional=True, description=
                        "report spliced alignment with score s only if <f>*s is "
                        "larger than next best spliced alignment "
                        "(default:1.000000)")
        self.add_option('hitstrategy', int, choices=[0, 1], optional=True,
                        default=1, description="report only best scoring hits "
                        "(=1) or all (=0) (default:1)")
        self.add_option('showalign', bool, optional=True, description=
                        "show alignments")
        self.add_option('prime5', str, optional=True, description=
                        "add 5' adapter (default:none)")
        self.add_option('prime3', str, optional=True, description=
                        "add 3' adapter (default:none)")
        self.add_option('clipacc', int, optional=True, description=
                        "clipping accuracy (default:70)")
        self.add_option('order', bool, optional=True, description=
                        "sorts the output by chromsome and position (might take "
                        "a while!)")
        self.add_option('maxinsertsize', int, optional=True, description=
                        "maximum size of the inserts (paired end) "
                        "(default:5000)")

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "256k")

    # self - macht class-funktion draus.
    # run_ids_connections_files - hash : run id -> n connections -> m files
    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['index', 'index2', 'filebins','bisulfite',
                   'minsize','brief','checkidx','briefcigar','threads',
                   'readgroupfile','readgroupid','differences','jump',
                   'nosuflinks','evalue','maxsplitevalue','maxinterval',
                   'maxout','splits','MEOP','nohead','extensionpenalty',
                   'dropoff','accuracy','minsplicecover','minfragscore','minfraglen',
                   'splicescorescale','hitstrategy','showalign',
                   'prime5','prime3','clipacc','order','maxinsertsize']
#        options = ['bisulfite', 'minsize', 'brief', 'differences',
#                   'jump', 'evalue', 'maxsplitevalue', 'maxinterval', 'splits',
#                   'SEGEMEHL', 'MEOP', 'nohead', 'extensionscore',
#                   'extensionpenalty', 'dropoff', 'accuracy', 'minsplicecover',
#                   'minfragscore', 'minfraglen', 'splicescorescale',
#                   'hitstrategy', 'showalign', 'prime5', 'prime3', 'clipacc',
#                   'polyA', 'autoclip', 'hardclip', 'order', 'maxinsertsize']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                input_paths = [ y for x in [fr_input, sr_input] \
                               for y in x if y != None ]

                # Do we have paired end data?
                is_paired_end = False if sr_input == [None] else True

                if len(fr_input) != 1 or fr_input == [None]:
                    logger.error("Expected single input file for first read.")
                    sys.exit(1)
                if is_paired_end and len(sr_input) != 1:
                    logger.error("Expected single input file for second read.")
                    sys.exit(1)

                if not os.path.isfile(self.get_option('index')):
                    logger.error(
                        "The path %s provided to option 'index' is not a file."
                        % self.get_option('index') )
                    sys.exit(1)

                if not os.path.isfile(self.get_option('index2')):
                    logger.error(
                        "The path %s provided to option 'index2' is not a file."
                        % self.get_option('index2') )
                    sys.exit(1)

                if not os.path.isfile(self.get_option('database')):
                    logger.error(
                        "The path %s provided to option 'database' is not a file."
                        % self.get_option('database'))
                    sys.exit(1)
                # SEGEMEHL can cope with gzipped files so we do not need to!!!
                #is_fr_gzipped = True if os.path.splitext(first_read_file[0])[1]\
                #                 in ['.gz', '.gzip'] else False

                #is_sr_gzipped = True if os.path.splitext(second_read_file[0])[1]\
                #                 in ['.gz', '.gzip'] else False

                # Segemehl is run in this exec group
                # Can segemehl handle multiple input files/fifos?
                with run.new_exec_group() as exec_group:
                    # 1. Create FIFO for genome
                    fifo_path_genome = run.add_temporary_file(
                        'segemehl-genome-fifo', designation = 'input')
                    mkfifo_genome = [self.get_tool('mkfifo'), fifo_path_genome]
                    exec_group.add_command(mkfifo_genome)
                    # 2. Create FIFO to write unmapped reads to
                    fifo_path_unmapped = run.add_temporary_file(
                        'segemehl-unmapped-fifo', designation = 'output')
                    mkfifo_unmapped = [
                        self.get_tool('mkfifo'),
                        fifo_path_unmapped
                    ]
                    exec_group.add_command(mkfifo_unmapped)
                    # 3. Read genome and output to FIFO
                    dd_genome = [self.get_tool('dd'),
                                 'bs=%s' % self.get_option('dd-blocksize'),
                                 'if=%s' % self.get_option('genome'),
                                 'of=%s' % fifo_path_genome]
                    exec_group.add_command(dd_genome)
                
                    with exec_group.add_pipeline() as segemehl_pipe:
                        # 4. Start segemehl
                        segemehl = [
                            self.get_tool('segemehl'),
                            '--database', fifo_path_genome,
                            '--nomatchfilename', fifo_path_unmapped,
                            '--threads', str(self.get_option('threads')),
                            '--query', fr_input[0]
                        ]
                        if is_paired_end:
                            segemehl.extend(['--mate', sr_input[0]])
                        segemehl.extend(option_list)
                        segemehl_pipe.add_command(
                            segemehl,
                            stderr_path = run.add_output_file(
                                'log',
                                '%s-segemehl-log.txt' % run_id,
                                input_paths)
                        )
                        # 4.1 command: Fix QNAMES in input SAM, if need be
                        if self.get_option('fix-qnames'):
                            fix_qnames = [
                                self.get_tool('fix_qnames'),
                                '--filetype', 'SAM'
                            ]
                            segemehl_pipe.add_command(fix_qnames)

                        # 5. Compress segemehl mapped reads
                        pigz_mapped_reads = [
                            self.get_tool('pigz'),
                            '--stdout',
                            '--blocksize', '4096', 
                            '--processes', '2'
                        ]

                        segemehl_pipe.add_command(
                            pigz_mapped_reads,
                            stdout_path = run.add_output_file(
                                'alignments',
                                '%s-segemehl-results.sam.gz' % run_id,
                                input_paths)
                        )

                    with exec_group.add_pipeline() as compress_unmapped_pipe:

                        # 6. Read unmapped reads from fifo
                        cat_unmapped_reads = [self.get_tool('cat'),
                                              fifo_path_unmapped]
                        compress_unmapped_pipe.add_command(cat_unmapped_reads)
                        
                        # 6.1 command: Fix QNAMES in input SAM, if need be
                        if self.get_option('fix-qnames'):
                            fix_qnames = [
                                self.get_tool('fix_qnames'),
                                '--filetype', 'FASTQ'
                            ]
                            compress_unmapped_pipe.add_command(fix_qnames)

                        # 7. Compress unmapped reads
                        pigz_unmapped_reads = [
                            self.get_tool('pigz'),
                            '--stdout',
                            '--blocksize', '4096',
                            '--processes', '1'
                        ]
                        compress_unmapped_pipe.add_command(
                            pigz_unmapped_reads,
                            stdout_path = run.add_output_file(
                                'unmapped',
                                '%s-segemehl-unmapped.fastq.gz' % run_id,
                                input_paths)
                        )
