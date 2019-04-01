import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class StringTie(AbstractStep):

    '''StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step
    to assemble and quantitate full-length transcripts representing multiple splice variants for
    each gene locus. Its input can include not only the alignments of raw reads used by other
    transcript assemblers, but also alignments longer sequences that have been assembled from those
    reads.In order to identify differentially expressed genes between experiments, StringTie's
    output can be processed by specialized software like Ballgown, Cuffdiff or other programs
    (DESeq2, edgeR, etc.).

    NOTE: This step implements that part of stringtie that assembles new transcripts. If you want
    stringtie to assemble transcripts from multiple input files please use step stringtie_merge!

    https://ccb.jhu.edu/software/stringtie/

    '''
    def __init__(self, pipeline):
        super(StringTie, self).__init__(pipeline)

        self.set_cores(6)

        # The BAM files
        self.add_connection('in/alignments',
                            constraints ={'min_files_per_run': 1, 'max_files_per_run': 1})
        # A .gtf file used as guide for assembly
        self.add_connection('in/features',
                            constraints = {'total_files': 1})

        self.add_connection('out/features')   # contains the assempled transcripts (GTF), -o
        self.add_connection('out/abundances') # -A <FILE.tab>
        self.add_connection('out/covered')    # -C <FILE.gtf>, requires -G!
#        self.add_connection('out/ballgown')   # -B, requires -G!
        # ballgown files
        # I skip this option if I don't know how many output files are returned here
        # this is maintained via uap anyway
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.require_tool('stringtie')
        self.require_tool('mkfifo')
        self.require_tool('dd')

        ## options for stringtie program
        # -G <FILE.gtf/gff>
        self.add_option('G', str, optional = True,
                        description = 'reference annotation to use for guiding the assembly process '
                        '(GTF/GFF3)')
        self.add_option('library_type', str, optional=False,
                        choices=['rf', 'fr'],
                        description='Assume stranded library fr-firststrand (--rf) or fr-secondstrand '
                        '(--fr).')
        # -l <LABEL>
        self.add_option('l', str, optional = True,
                        description = 'name prefix for output transcripts (default: STRG)')
        # -f <0.1-1.0>
        self.add_option('f', float, optional = True,
                        description = 'minimum isoform fraction (default: 0.1)')
        # -m <INT>
        self.add_option('m', int, optional = True,
                        description = 'minimum assembled transcript length (default: 200)')
        # -a <INT>
        self.add_option('a', int, optional = True,
                        description = 'minimum anchor length for junctions (default: 10)')
        # -j <FLOAT>
        self.add_option('j', float, optional = True,
                        description = 'minimum junction coverage (default: 1)')
        # -t
        self.add_option('t', bool, optional = True,
                        description = 'disable trimming of predicted transcripts based on coverage '
                        '(default: coverage trimming is enabled)')
        # -c <FLOAT>
        self.add_option('c', float, optional = True,
                        description = 'minimum reads per bp coverage to consider for transcript '
                        'assembly (default: 2.5)')
        # -v
        self.add_option('v', bool, optional=True,
                        description='verbose (log bundle processing details)')
        # -g <INT>
        self.add_option('g', int, optional = True,
                        description = 'gap between read mappings triggering a new bundle (default: 50)')
        # -C <FILE.gtf> => out/covered, requires -G -> might be an empty file if !-G 
	# this file is generated in any case 
#        self.add_option('covered-references', bool, optional = True,
#                        description = 'Write reference transcripts that are covered by reads to an '
#                        'output .gtf file. The file will be empty if this option is not set. Default: False')
        # -M <0.0-1.0>
        self.add_option('M', float, optional = True,
                        description = 'fraction of bundle allowed to be covered by multi-hit reads '
                        '(default:0.95)')
        # -p <INT>
        self.add_option('p', int, optional=True,
                        description='number of threads (CPUs) to use (default: 1)')
        # OLLI -> I specified 4 in a test run and only 2 were used!
        # Q: is it only capable of using 2 CPUs?

        # -A <FILE.tab> => out/abundances
	# this file is generated in any case via uap
#        self.add_option('abundances', bool, optional = True,
#                        description = 'Print gene abundance estimation to an output file. The file will '
#                        'be empty if this option is not set. Default: False')
        # ballgown files
        self.add_option('B', bool, optional = True,
                        description = 'Enable the ouput of Ballgown input table files (.ctab). '
                        'containing coverage data for the reference transcripts given with the '
                        '-G option. (See the Ballgown documentation for a description of these '
                        'files.) With this option StringTie can be used as a direct replacement '
                        'of the tablemaker program included with the Ballgown distribution. ')
        # -B, -b this is covered with -b option that specifies the exact path for this file
        # If the option -o is given as a full path to the output transcript file, StringTie
        # will write the *.ctab files in the same directory as the output GTF.
        # -b option is not available for uap, since it decides on its own where the output will be stored.

        # -e
        self.add_option('e', bool, optional = True, default = False,
                        description = 'only estimate the abundance of given reference transcripts '
                        '(requires -G)')
        # -x <seqid_list>
        self.add_option('x', str, optional = True,
                        description = 'Ignore all read alignments (and thus do not attempt to '
                        'perform transcript assembly) on the specified reference sequences. '
                        'Parameter <seqid_list> can be a single reference sequence name (e.g. '
                        '-x chrM) or a comma-delimited list of sequence names (e.g. -x '
                        '"chrM,chrX,chrY"). This can speed up StringTie especially in the case '
                        'of excluding the mitochondrial genome, whose genes may have very high '
                        'coverage in some cases, even though they may be of no interest for a '
                        'particular RNA-Seq analysis. The reference sequence names are case '
                        'sensitive, they must match identically the names of chromosomes/contigs '
                        'of the target genome against which the RNA-Seq reads were aligned in '
                        'the first place.')

        # OLLI --merge, needs an additional uap step since it has different input/output conncections.

        self.add_option('merged-id', str, optional=True, default="mergeSBam",
                        description="The run-id that has been used in the samtools_merge step before")

        # [Options for 'dd':]
	self.add_option('fifo', bool, optional = True, default = False,
			description='Enable the FIFO functionality for splitting large input files.')	
        self.add_option('dd-blocksize', str, optional = True, default = "2M",
			description='Provide the blocksize for dd tool.')

    def runs(self, run_ids_connections_files):

        # Compile the list of options
        options=['l','f','m','a','j','t','c','v','g','M','p', 'B','e','x', 'G']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append( '-%s' % option )
                option_list.append( str(self.get_option(option)) )

        # library-type
        if self.is_option_set_in_config('library_type'):
            option_list.append('--%s' % self.get_option('library_type'))

        merged_id = self.get_option('merged-id')

        features_path = str

        for run_id in run_ids_connections_files.keys():
            try:
                features_path = run_ids_connections_files[run_id]['in/features'][0]
            except KeyError:
                if self.is_option_set_in_config('G'):
                    features_path = self.get_option('G')
                #else:
                #    logger.error("No feature file could be found for '%s'" % run_id)
                #    sys.exit(1)

        #if not features_path:
        #    option_list.append('-G')
        #    option_list.append(features_path)

        for run_id in run_ids_connections_files.keys():

            try:
                input_paths = run_ids_connections_files[run_id]['in/alignments']
            except KeyError:
                continue # this happens if merged annotation from previous step is used as reference
                         # guide

            with self.declare_run(run_id) as run:
                outfile = run.add_output_file('features',
                                              '%s-transcripts.gtf' % run_id,
                                              input_paths)
                abundfile = run.add_output_file('abundances',
                                                '%s-abundances.gtf' % run_id,
                                                input_paths)
                covfile = run.add_output_file('covered',
                                              '%s-coveredRefs.gtf' % run_id,
                                              input_paths)
                stdout = run.add_output_file('log_stdout',
                                             '%s-stringtie.stdout' % run_id,
                                             input_paths)
                stderr = run.add_output_file('log_stderr',
                                             '%s-stringtie.stderr' % run_id,
                                             input_paths)


                with run.new_exec_group() as exec_group:

                    if self.get_option('fifo'):
			# 1. create FIFO for BAM file
	                fifo_path_bam = run.add_temporary_file('bam_path_fifo',
        	                                               designation = 'input')
                    	mkfifo_bam = [self.get_tool('mkfifo'), fifo_path_bam]
                    	exec_group.add_command(mkfifo_bam)

                    	# 2. read BAM and output to FIFO
                    	dd_bam = [self.get_tool('dd'),
                                  'bs=%s' % self.get_option('dd-blocksize'),
                              	  'if=%s' % input_paths[0],
                                  'of=%s' % fifo_path_bam]
                    	exec_group.add_command(dd_bam)

                   	# 3. initialize the stringtie command on FIFO
                    	stringtie = [self.get_tool('stringtie'), fifo_path_bam]
		    else:
			# 1. initialize the stringtie command on input BAM file
                    	stringtie = [self.get_tool('stringtie'), input_paths[0]]

                    stringtie.extend(option_list)

 	            stringtie.extend(['-C', covfile])

                    stringtie.extend(['-A', abundfile])

                    stringtie.extend(['-o', outfile])

                    exec_group.add_command(stringtie,
                                           stdout_path = stdout,
                                           stderr_path = stderr)
