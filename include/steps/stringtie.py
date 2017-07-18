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

    https://ccb.jhu.edu/software/stringtie/

    '''
    def __init__(self, pipeline):
        super(StringTie, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/features')   # contains the assempled transcripts (GTF), -o
        self.add_connection('out/abundances') # -A <FILE.tab>
        self.add_connection('out/covered')    # -C <FILE.gtf>, requires -G!
#        self.add_connection('out/ballgown')   # -B, requires -G!
        # ballgown files
        # I skip this option if I don't know how many output files are returned here
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')
        
        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('stringtie')

        ## options for stringtie program
        # -p <INT>
        self.add_option('p', int, optional=True, default = 1,
                        description='number of threads used during analysis')
        # -> I specified 4 in a test run and only 2 were used!
        # Q: is it only capable of using 2 CPUs?

        # -G <FILE.gtf/gff>
        self.add_option('G', str, optional = True,
                        description = 'Use the reference annotation file (in GTF or GFF3 format) '
                        'to guide the assembly process. The output will include expressed reference '
                        'transcripts as well as any novel transcripts that are assembled. This '
                        'option is required by options -B, -b, -e, -C (see below). ')
        # --rf/--fr
        self.add_option('library-type', str, optional=False,choices=['fr-firststrand', 'fr-secondstrand'], 
                        description='Assumes a stranded library. Allowed values are "fr-firststrand '
                        'and "fr-secondstrand".')
        # -l <LABEL>
        self.add_option('l', str, optional = True, default = "STRG",
                        description = 'Sets <LABEL> as the prefix for the name of the output '
                        'transcripts.') 
        # -f <0.1-1.0>
        self.add_option('f', float, optional = True, default = 0.1,
                        description = 'Sets the minimum isoform abundance of the predicted '
                        'transcripts as a fraction of the most abundant transcript assembled at '
                        'a given locus. Lower abundance transcripts are often artifacts of '
                        'incompletely spliced precursors of processed transcripts. Values must'
                        'be in the interval <0.1-1.0>.')
        # -m <INT>
        self.add_option('m', int, optional = True, default = 200,
                        description = 'Sets the minimum length allowed for the predicted '
                        'transcripts.')

        # -A <FILE.tab> => out/abundances
        self.add_option('abundances', bool, optional = True, default = False,
                        description = 'Print a table with the gene abundances.')

        # -C <FILE.gtf> => out/covered, requires -G -> might be an empty file if !-G 
        self.add_option('covered-references', bool, optional = True, default = False,
                        description = 'Print a gtf file of all the covered reference transcripts.')

        # -a <INT>
        self.add_option('a', int, optional = True, default = 10, 
                        description = 'Junctions that do not have spliced reads that align '
                        'across them with at least this amount of bases on both sides are '
                        'filtered out.')
        # -j <FLOAT>
        self.add_option('j', float, optional = True, default = 1.0,
                        description = 'There should be at least this many spliced reads that '
                        'align across a junction (i.e. junction coverage). This number can be '
                        'fractional, since some reads align in more than one place. A read that '
                        'aligns in n places will contribute 1/n to the junction coverage.')
        # -t
        self.add_option('t', bool, optional = True, default = False,
                        description = 'This parameter disables trimming at the ends of the '
                        'assembled transcripts. By default StringTie adjusts the predicted '
                        'transcripts start and/or stop coordinates based on sudden drops in '
                        'coverage of the assembled transcript. ') 
        # -c <FLOAT>
        self.add_option('c', float, optional = True, default = 2.5,
                        description = 'Sets the minimum read coverage allowed for the predicted '
                        'transcripts. A transcript with a lower coverage than this value is not '
                        'shown in the output.')
        # -g <INT>
        self.add_option('g', int, optional = True, default = 50,
                        description = 'Minimum locus gap separation value (in bp). Reads that '
                        'are mapped closer than this distance are merged together in the same '
                        'processing bundle.')

        # ballgown files
        self.add_option('ballgown', bool, optional = True, default = False,
                        description = 'Enable the ouput of Ballgown input table files (.ctab). '
                        'containing coverage data for the reference transcripts given with the '
                        '-G option. (See the Ballgown documentation for a description of these '
                        'files.) With this option StringTie can be used as a direct replacement '
                        'of the tablemaker program included with the Ballgown distribution. ')
        # -B, -b this is covered with -b option that specifies the exact path for this file
        # If the option -o is given as a full path to the output transcript file, StringTie 
        # will write the *.ctab files in the same directory as the output GTF.

        # -e
        self.add_option('e', bool, optional = True, default = False,
                        description = 'Limits the processing of read alignments to only estimate '
                        'and output the assembled transcripts matching the reference transcripts '
                        'given with the -G option (requires -G, recommended for -B/-b). With this '
                        'option, read bundles with no reference transcripts will be entirely '
                        'skipped, which may provide a considerable speed boost when the given set '
                        'of reference transcripts is limited to a set of target genes, for example.')
        
        # -M <0.0-1.0>
        self.add_option('M', float, optional = True, default = 0.95,
                        description = 'Sets the maximum fraction of muliple-location-mapped '
                        'reads that are allowed to be present at a given locus.')

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

        # --merge, needs an additional uap step since it has different input/output conncections.

    def runs(self, run_ids_connections_files):
        
        # Compile the list of options
        options=['p', 'G', 'l', 'f', 'm', 'a,', 'j', 't', 'c', 'g', 'e', 'M', 'x']

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
        if self.get_option('library-type') == 'fr-firststrand':
            option_list.append('--rf')
        elif self.get_option('library-type') == 'fr-secondstrand':
            option_list.append('--fr')
        else:
            raise StandardError('Unexpected value for option "library-type" in step "stringtie". '
                                'Accepted values are: "fr-firststrand" or "fr-secondstrand".')

        for run_id in run_ids_connections_files.keys():

             with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/alignments']
#                temp_dir = run.add_temporary_directory('cufflinks-out')
                
                # check, if only a single input file is provided
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file., but got this %s" % input_paths)

                outfile = run.add_output_file('features', '%s-transcripts.gtf' % run_id, input_paths)
                abundfile = run.add_output_file('abundances', '%s-abundances.gtf' % run_id, input_paths)
                covfile = run.add_output_file('covered', '%s-coveredRefs.gtf' % run_id, input_paths)
                stdout = run.add_output_file('log_stdout', '%s-stringtie.stdout' % run_id, input_paths)
                stderr = run.add_output_file('log_stderr', '%s-stringtie.stderr' % run_id, input_paths)

                stringtie = [self.get_tool('stringtie'), input_paths[0]]
                stringtie.extend(option_list)

                if self.get_option('abundances'):
                    stringtie.extend(['-A', abundfile])

                if self.get_option('covered-references'):
                    stringtie.extend(['-C', covfile])

                stringtie.extend(['-o', outfile])

                with run.new_exec_group() as exec_group:
                    exec_group.add_command(stringtie,
                                           stderr_path = stderr,
                                           stdout_path = stdout)
