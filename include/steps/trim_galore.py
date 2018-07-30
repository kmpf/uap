import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')

class TrimGalore(AbstractStep):
    '''
    A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to
    FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation
    Bisufite-Seq) libraries.

    https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

    Note for RRBS using the NuGEN Ovation RRBS System 1-16 kit:

    Owing to the fact that the NuGEN Ovation kit attaches a varying number of nucleotides (0-3) after each MspI
    site Trim Galore should be run WITHOUT the option --rrbs. This trimming is accomplished in a subsequent 
    diversity trimming step afterwards (see their manual).

    Note for RRBS using MseI:

    If your DNA material was digested with MseI (recognition motif: TTAA) instead of MspI it is NOT necessary
    to specify --rrbs or --non_directional since virtually all reads should start with the sequence
    'TAA', and this holds true for both directional and non-directional libraries. As the end-repair of 'TAA'
    restricted sites does not involve any cytosines it does not need to be treated especially. Instead, simply
    run Trim Galore! in the standard (i.e. non-RRBS) mode.

    '''

    def __init__(self, pipeline):

        super(TrimGalore, self).__init__(pipeline)

        self.set_cores(4) # consider to adjust this

# Writing report to 'MAIT_0h_A_R1.fastq.gz_trimming_report.txt'
# Writing final adapter and quality trimmed output to MAIT_0h_A_R1_trimmed.fq.gz
# Writing to stdout and stderr

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        self.add_connection('out/first_read') # <run_id>_R1_val_1.fq.gz
        self.add_connection('out/first_read_report') #_trimming_report.txt
        self.add_connection('out/first_read_fastqc_zip') # <run_id>_R1_val_1.fq.gz
        self.add_connection('out/first_read_fastqc_html') # <run_id>_R1_val_1.fq.html

        self.add_connection('out/second_read')
        self.add_connection('out/second_read_report') #_trimming_report.txt
        self.add_connection('out/second_read_fastqc_zip') # <run_id>_R1_val_1.fq.gz
        self.add_connection('out/second_read_fastqc_html') # <run_id>_R1_val_1.fq.html

        self.add_connection('out/stdout')
        self.add_connection('out/stderr')

        self.require_tool('trim_galore')
        self.require_tool('cutadapt')

# - copy help message into python file
# - rm newlines from descriptions via Alt+Ctrl+Shift+5 RET C-q C-j RET RET
# --> this replaces nl with 25 whitespace
# - rm this whitespace via M-x delete-horizontal-space
# rm all blank lines
# - add self.add_option(' before each line
# - replace types to python types (bool for all non-type options!)
# - add description=' before description
# - add ') to each line end
# - insert optional=True, before the word description
# - replace whitespace between type and optional=True via Alt+Ctrl+Shift+5 RET \s-+optional= RET  optional=
# - insert lb before description
# - remove all trailing whitespace via M-x delete-trailing-whitespace

        self.add_option('quality', int, optional=True,
                        description='Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming'
                        'will be performed first, and adaptertrimming is carried in a second round. Other files are quality'
                        'and adaptertrimmed in a single pass. The algorithm is the same as the one used by BWA(Subtract INT'
                        'from all qualities; compute partial sums from all indicesto the end of the sequence; cut sequence at'
                        'the index at which the sum isminimal). Default Phred score: 20.')
        self.add_option('phred33', bool, optional=True,
                        description='Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding)'
                        'for quality trimming. Default: ON.')
        self.add_option('phred64', bool, optional=True,
                        description='Instructs Cutadapt to use ASCII+64 quality scores as Phred scores(Illumina 1.5 encoding) for'
                        'quality trimming.')
        self.add_option('fastqc', bool, optional=True,
                        description='Run FastQC in the default mode on the FastQ file once trimming is complete.')
        self.add_option('fastqc_args', str, optional=True,
                        description='Passes extra arguments to FastQC. If more than one argument is to be passedto FastQC they must be'
                        'in the form "arg1 arg2 etc.". An example would be:--fastqc_args "--nogroup --outdir /home/". Passing'
                        'extra arguments willautomatically invoke FastQC, so --fastqc does not have to be'
                        'specifiedseparately.')
        self.add_option('adapter', str, optional=True,
                        description='Adapter sequence to be trimmed. If not specified explicitly, Trim Galore willtry to auto-detect'
                        'whether the Illumina universal, Nextera transposase or Illuminasmall RNA adapter sequence was used.'
                        'Also see \'--illumina\', \'--nextera\' and \'--small_rna\'. If no adapter can be detected within the'
                        'first 1 million sequencesof the first file specified Trim Galore defaults to \'--illumina\'.')
        self.add_option('adapter2', str, optional=True,
                        description='Optional adapter sequence to be trimmed off read 2 of paired-end files. Thisoption requires'
                        '\'--paired\' to be specified as well. If the libraries to be trimmedare smallRNA then a2 will be set'
                        'to the Illumina small RNA 5\' adapter automatically(GATCGTCGGACT).')
        self.add_option('illumina', bool, optional=True,
                        description='Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter'
                        '\'AGATCGGAAGAGC\' instead of the default auto-detection of adapter sequence.')
        self.add_option('nextera', bool, optional=True,
                        description='Adapter sequence to be trimmed is the first 12bp of the Nextera adapter \'CTGTCTCTTATA\' instead of'
                        'the default auto-detection of adapter sequence.')
        self.add_option('small_rna',bool, optional=True,
                        description='Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3\' Adapter'
                        '\'TGGAATTCTCGG\' instead of the default auto-detection of adapter sequence. Selectingto trim'
                        'smallRNA adapters will also lower the --length value to 18bp. If the smallRNAlibraries are'
                        'paired-end then a2 will be set to the Illumina small RNA 5\' adapterautomatically (GATCGTCGGACT)'
                        'unless -a 2 had been defined explicitly.')
        self.add_option('max_length', int, optional=True,
                        description='Discard reads that are longer than <INT> bp after trimming. This is only advised for smallRNA'
                        'sequencing to remove non-small RNA sequences.')
        self.add_option('stringency', int, optional=True,
                        description='Overlap with adapter sequence required to trim a sequence. Defaults to avery stringent setting of'
                        '1, i.e. even a single bp of overlapping sequencewill be trimmed off from the 3\' end of any read.')
        self.add_option('e', float, optional=True,
                        description='Maximum allowed error rate (no. of errors divided by the length of the matchingregion) (default:'
                        '0.1)')
        self.add_option('gzip', bool, optional=True,
                        description='Compress the output file with GZIP. If the input files are GZIP-compressedthe output files will'
                        'automatically be GZIP compressed as well. As of v0.2.8 the compression will take place on the fly.')
        self.add_option('dont_gzip', bool, optional=True,
                        description='Output files won\'t be compressed with GZIP. This option overrides --gzip.')
        self.add_option('length', int, optional=True,
                        description='Discard reads that became shorter than length INT because of eitherquality or adapter trimming. A'
                        'value of \'0\' effectively disablesthis behaviour. Default: 20 bp.For paired-end files, both reads'
                        'of a read-pair need to be longer than<INT> bp to be printed out to validated paired-end files (see'
                        'option --paired).If only one read became too short there is the possibility of keeping suchunpaired'
                        'single-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.')
        self.add_option('max_n', int, optional=True,
                        description='The total number of Ns (as integer) a read may contain before it will be removed altogether.In a'
                        'paired-end setting, either read exceeding this limit will result in the entirepair being removed'
                        'from the trimmed output files.')
        self.add_option('trim-n', bool, optional=True,
                        description='Removes Ns from either side of the read. This option does currently not work in RRBS mode.')
        # output dir is defined by uap (a la cuffcompare), no need to add this option here
#        self.add_option('-o/--output_dir', str, optional=True,
#                        description='If specified all output will be written to this directory instead of the currentdirectory.')
        # uap will always expect a report file, thus, this option is not available for uap users
#        self.add_option('no_report_file', bool, optional=True,
#                        description='If specified no report file will be generated.')
        # uap will always report warnings to stdout and stderr, thus, this option is not available for uap users
#        self.add_option('suppress_warn',bool, optional=True,
#                        description='If specified any output to STDOUT or STDERR will be suppressed.')
        self.add_option('clip_R1', int, optional=True,
                        description='Instructs Trim Galore to remove <int> bp from the 5\' end of read 1 (or single-endreads). This may'
                        'be useful if the qualities were very poor, or if there is somesort of unwanted bias at the 5\' end.'
                        'Default: OFF.')
        self.add_option('clip_R2', int, optional=True,
                        description='Instructs Trim Galore to remove <int> bp from the 5\' end of read 2 (paired-end readsonly). This'
                        'may be useful if the qualities were very poor, or if there is some sortof unwanted bias at the 5\''
                        'end. For paired-end BS-Seq, it is recommended to removethe first few bp because the end-repair'
                        'reaction may introduce a bias towards lowmethylation. Please refer to the M-bias plot section in the'
                        'Bismark User Guide forsome examples. Default: OFF.')
        self.add_option('three_prime_clip_R1', int, optional=True,
                        description='Instructs Trim Galore to remove <int> bp from the 3\' end of read 1 (or single-endreads) AFTER'
                        'adapter/quality trimming has been performed. This may remove some unwantedbias from the 3\' end that'
                        'is not directly related to adapter sequence or basecall quality.Default: OFF.')
        self.add_option('three_prime_clip_R2', int, optional=True,
                        description='Instructs Trim Galore to remove <int> bp from the 3\' end of read 2 AFTERadapter/quality trimming'
                        'has been performed. This may remove some unwanted bias fromthe 3\' end that is not directly related'
                        'to adapter sequence or basecall quality.Default: OFF.')
        # uap requires the tool cutadapt, so it will be in the path and this option is not necessary for uap users
#        self.add_option('path_to_cutadapt', str, optional=True,
#                        description='You may use this option to specify a path to the Cutadapt executable,e.g.'
#                        '/my/home/cutadapt-1.7.1/bin/cutadapt. Else it is assumed that Cutadapt is in the PATH.')

        # RRBS-specific options (MspI digested material):
        self.add_option('rrbs', bool, optional=True,
                        description='Specifies that the input file was an MspI digested RRBS sample (recognitionsite: CCGG). Single-end'
                        'or Read 1 sequences (paired-end) which were adapter-trimmedwill have a further 2 bp removed from'
                        'their 3\' end. Sequences which were merelytrimmed because of poor quality will not be shortened'
                        'further. Read 2 of paired-endlibraries will in addition have the first 2 bp removed from the 5\' end'
                        '(by setting \'--clip_r2 2\'). This is to avoid using artificial methylation calls from the'
                        'filled-incytosine positions close to the 3\' MspI site in sequenced fragments. This option is not'
                        'recommended for users of the NuGEN ovation RRBS System 1-16kit (see below).')
        self.add_option('non_directional', bool, optional=True,
                        description='Selecting this option for non-directional RRBS libraries will screenquality-trimmed sequences for'
                        '\'CAA\' or \'CGA\' at the start of the readand, if found, removes the first two basepairs. Like with'
                        'the option \'--rrbs\' this avoids using cytosine positions that were filled-induring the end-repair'
                        'step. \'--non_directional\' requires \'--rrbs\' tobe specified as well. Note that this option does'
                        'not set \'--clip_r2 2\' inpaired-end mode.')
        self.add_option('keep', bool, optional=True,
                        description='Keep the quality trimmed intermediate file. Default: off, which meansthe temporary file is being'
                        'deleted after adapter trimming. Only hasan effect for RRBS samples since other FastQ files are not'
                        'trimmedfor poor qualities separately.')

        # Paired-end specific options:
        self.add_option('paired', bool, optional=True,
                        description='This option performs length trimming of quality/adapter/RRBS trimmed reads forpaired-end files. To'
                        'pass the validation test, both sequences of a sequence pairare required to have a certain minimum'
                        'length which is governed by the option--length (see above). If only one read passes this length'
                        'threshold theother read can be rescued (see option --retain_unpaired). Using this option letsyou'
                        'discard too short read pairs without disturbing the sequence-by-sequence orderof FastQ files which'
                        'is required by many aligners.Trim Galore! expects paired-end files to be supplied in a pairwise'
                        'fashion, e.g.file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... .')
        self.add_option('trim1', bool, optional=True,
                        description='Trims 1 bp off every read from its 3\' end. This may be needed for FastQ files thatare to be'
                        'aligned as paired-end data with Bowtie. This is because Bowtie (1) regardsalignments like this: R1'
                        '---------------------------> or this: -----------------------> R1 R2 <---------------------------'
                        '<----------------- R2as invalid (whenever a start/end coordinate is contained within the other'
                        'read).NOTE: If you are planning to use Bowtie2, BWA etc. you don\'t need to specify this option.')
        self.add_option('retain_unpaired', bool, optional=True,
                        description='If only one of the two paired-end reads became too short, the longerread will be written to either'
                        '\'.unpaired_1.fq\' or \'.unpaired_2.fq\' output files. The length cutoff for unpaired single-end'
                        'reads isgoverned by the parameters -r1/--length_1 and -r2/--length_2. Default: OFF.')
        self.add_option('length_1', int, optional=True,
                        description='Unpaired single-end read length cutoff needed for read 1 to be written to \'.unpaired_1.fq \''
                        'output file. These reads may be mapped in single-end mode.Default: 35 bp.')
        self.add_option('length_2', int, optional=True,
                        description='Unpaired single-end read length cutoff needed for read 2 to be written to \'.unpaired_2.fq \''
                        'output file. These reads may be mapped in single-end mode.Default: 35 bp.')

        # any other options for uap behaviour?
        # [Options for 'dd' and 'pigz':]
        self.add_option('dd-blocksize', str, optional = True, default = "2M")
        self.add_option('pigz-blocksize', str, optional = True, default = "2048")

    def runs(self, run_ids_connections_files):

        options = ['quality','phred33','phred64','fastqc','adapter','adapter2',
                   'illumina','nextera','small_rna','max_length','stringency', 'gzip','dont_gzip',
                   'length','max_n','trim-n', 'clip_R1','clip_R2','three_prime_clip_R1',
                   'three_prime_clip_R2', 'rrbs','non_directional','keep',
                   'paired','trim1','retain_unpaired','length_1','length_2']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append( '--%s' % option )
                option_list.append( str(self.get_option(option)) )

        # this is the only option that only a '-' instead of '--' in the call
        option = 'e'
        if self.is_option_set_in_config(option):
            option_list.append('-%s' % option)
            option_list.append(str(self.get_option(option)))

        # fastqc_args needs extra quote characters
        option = 'fastqc_args'
        if self.is_option_set_in_config(option):
            option_list.append('--%s' % option)
            option_list.append('\"%s\"' % str(self.get_option(option)))

        # set possible read types
        read_types = {'first_read': 'R1', 'second_read': 'R2'}
        paired_end_info = dict()

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:

                # the temporary output directory
                outdir = run.get_output_directory_du_jour_placeholder()
                # this is the prefix for the trim_galore cmd option:
                # -o
                prefixTG = '%s' % outdir

                input_paths = run_ids_connections_files[run_id]

                if input_paths['in/first_read'][0].endswith('.gz'):
                    run.add_output_file('first_read',
                                        '%s_R1_val_1.fq.gz' % run_id,
                                        input_paths['in/first_read'])
                else:
                    run.add_output_file('first_read',
                                        '%s_R1_val_1.fq' % run_id,
                                        input_paths['in/first_read'])

                run.add_output_file('first_read_fastqc_zip',
                                    '%s_R1_val_1_fastqc.zip' % run_id,
                                    input_paths['in/first_read'])
                run.add_output_file('first_read_fastqc_html',
                                    '%s_R1_val_1_fastqc.html' % run_id,
                                    input_paths['in/first_read'])
                run.add_output_file('first_read_report',
                                    '%s_trimming_report.txt' % os.path.basename(input_paths['in/first_read'][0]),
                                    input_paths['in/first_read'])

                if input_paths['in/second_read'] :
                    if input_paths['in/second_read'][0].endswith('.gz'):
                        run.add_output_file('second_read',
                                            '%s_R2_val_2.fq.gz' % run_id,
                                            input_paths['in/second_read'])
                    else:
                        run.add_output_file('second_read',
                                            '%s_R2_val_2.fq' % run_id,
                                            input_paths['in/second_read'])

                    run.add_output_file('second_read_fastqc_zip',
                                        '%s_R2_val_2_fastqc.zip' % run_id,
                                        input_paths['in/second_read'])
                    run.add_output_file('second_read_fastqc_html',
                                        '%s_R2_val_2_fastqc.html' % run_id,
                                        input_paths['in/second_read'])
                    run.add_output_file('second_read_report',
                                        '%s_trimming_report.txt' % os.path.basename(input_paths['in/second_read'][0]),
                                        input_paths['in/second_read'])
                else:
                    run.add_empty_output_connection('second_read')
                    run.add_empty_output_connection('second_read_report')
                    run.add_empty_output_connection('second_read_fastqc_zip')
                    run.add_empty_output_connection('second_read_fastqc_html')

                run.add_output_file('first_read_report',
                                    '%s_R1_trimming_report.txt' % run_id,
                                    input_paths['in/first_read'])

                stdout = run.add_output_file('stdout',
                                             '%s_trimgalore.stdout' % run_id,
                                             input_paths['in/first_read'])
                stderr = run.add_output_file('stderr',
                                             '%s_trimgalore.stderr' % run_id,
                                             input_paths['in/first_read'])

                # set '--paired' option if not already been done by usr in cfg
                if input_paths['in/second_read']:
                    if not self.is_option_set_in_config('paired'):
                        option_list.append('--paired')



                tg = [self.get_tool('trim_galore'),
                              '--output_dir', outdir]
                tg.extend(option_list)
                tg.append(input_paths['in/first_read'][0])
                if input_paths['in/second_read']:
                    tg.append(input_paths['in/second_read'][0])

                with run.new_exec_group() as tg_exec_group:
                    tg_exec_group.add_command(tg,
                                              stderr_path = stderr,
                                              stdout_path = stdout)

