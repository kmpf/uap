from uaperrors import StepError
import sys
import os
import re
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class BwaBacktrack(AbstractStep):
    '''
    bwa-backtrack is the bwa algorithm designed for Illumina sequence reads up
    to 100bp. The computation of the alignments is done by running 'bwa aln'
    first, to align the reads, followed by running 'bwa samse' or 'bwa sampe'
    afterwards to generate the final SAM output.

    http://bio-bwa.sourceforge.net/

    typical command line for single-end data::

        bwa aln <bwa-index> <first-read.fastq> > <first-read.sai>
        bwa samse <bwa-index> <first-read.sai> <first-read.fastq> > <sam-output>

    typical command line for paired-end data::

        bwa aln <bwa-index> <first-read.fastq> > <first-read.sai>
        bwa aln <bwa-index> <second-read.fastq> > <second-read.sai>
        bwa sampe <bwa-index> <first-read.sai> <second-read.sai> \
                  <first-read.fastq> <second-read.fastq> > <sam-output>

    '''

    def __init__(self, pipeline):
        super(BwaBacktrack, self).__init__(pipeline)
        self.set_cores(8)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')

        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')
        # Step was tested for bwa release 0.7.15-r1140
        self.require_tool('bwa')

        # Options for the programs bwa aln/samse/sampe
        self.add_option('index', str, optional=False,
                        description="Path to BWA index")
        # [Options for 'bwa aln':]
        self.add_option(
            'aln-n',
            float,
            optional=True,
            description="Maximum edit distance if the value is "
            "INT, or the fraction of missing alignments given 2% "
            "uniform base error rate if FLOAT. In the latter case, "
            "the maximum edit distance is automatically chosen for "
            "different read lengths. [0.04]")
        self.add_option('aln-o', int, optional=True,
                        description="Maximum number of gap opens [1]")
        self.add_option('aln-e', int, optional=True,
                        description="Maximum number of gap extensions, -1 for "
                        "k-difference mode (disallowing long gaps) [-1]")
        self.add_option('aln-d', int, optional=True,
                        description="Disallow a long deletion within INT bp "
                        "towards the 3'-end [16]")
        self.add_option('aln-i', int, optional=True,
                        description="Disallow an indel within INT bp towards "
                        "the ends [5]")
        self.add_option(
            'aln-l',
            int,
            optional=True,
            description="Take the first INT subsequence as seed. "
            "If INT is larger than the query sequence, seeding will "
            "be disabled. For long reads, this option is typically "
            "ranged from 25 to 35 for '-k 2'. [inf]")
        self.add_option('aln-k', int, optional=True,
                        description="Maximum edit distance in the seed [2]")
        self.add_option('aln-t', int, optional=True, default=1,
                        description="Number of threads (multi-threading mode) "
                        "[1]")
        self.add_option('aln-M', int, optional=True,
                        description="Mismatch penalty. BWA will not search "
                        "for suboptimal hits with a score lower than "
                        "(bestScore-misMsc). [3]")
        self.add_option('aln-O', int, optional=True,
                        description="Gap open penalty [11]")
        self.add_option('aln-E', int, optional=True,
                        description="Gap extension penalty [4]")
        self.add_option('aln-R', int, optional=True,
                        description="Proceed with suboptimal alignments if "
                        "there are no more than INT equally best hits. This "
                        "option only affects paired-end mapping. Increasing "
                        "this threshold helps to improve the pairing accuracy "
                        "at the cost of speed, especially for short reads "
                        "(~32bp).")
        self.add_option('aln-c', bool, optional=True,
                        description="Reverse query but not complement it, "
                        "which is required for alignment in the color space. "
                        "(Disabled since 0.6.x)")
        self.add_option('aln-N', bool, optional=True,
                        description="Disable iterative search. All hits with "
                        "no more than maxDiff differences will be found. This "
                        "mode is much slower than the default.")
        self.add_option('aln-q', int, optional=True,
                        description="Parameter for read trimming. BWA trims a "
                        r"read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if "
                        "q_l<INT where l is the original read length. [0]")
        self.add_option('aln-I', bool, optional=True,
                        description="The input is in the Illumina 1.3+ read "
                        "format (quality equals ASCII-64).")
        self.add_option('aln-B', int, optional=True,
                        description="Length of barcode starting from the "
                        "5'-end. When INT is positive, the barcode of each read "
                        "will be trimmed before mapping and will be written at "
                        "the BC SAM tag. For paired-end reads, the barcode from "
                        "both ends are concatenated. [0]")
        self.add_option(
            'aln-b',
            bool,
            optional=True,
            description="Specify the input read sequence file is "
            "the BAM format. For paired-end data, two ends in a "
            "pair must be grouped together and options aln-1 or "
            "aln-2 are usually applied to specify which end should "
            "be mapped. Typical command lines for mapping pair-end "
            "data in the BAM format are:\n"
            "     bwa aln ref.fa -b1 reads.bam > 1.sai\n"
            "     bwa aln ref.fa -b2 reads.bam > 2.sai \n"
            "     bwa sampe ref.fa 1.sai 2.sai reads.bam reads.bam > "
            "aln.sam\n")
        self.add_option('aln-0', bool, optional=True,
                        description="When aln-b is specified, only use single-"
                        "end reads in mapping.")
        self.add_option(
            'aln-1',
            bool,
            optional=True,
            description="When aln-b is specified, only use the "
            "first read in a read pair in mapping (skip single-end "
            "reads and the second reads).")
        self.add_option('aln-2', bool, optional=True,
                        description="When aln-b is specified, only use the "
                        "second read in a read pair in mapping.")

        # [Options for 'bwa samse':]
        self.add_option(
            'samse-n',
            int,
            optional=True,
            description="Maximum number of alignments to output "
            "in the XA tag for reads paired properly. If a read "
            "has more than INT hits, the XA tag will not be written."
            " [3]")
        self.add_option('samse-r', str, optional=True,
                        description="Specify the read group in a format like "
                        "'@RG\tID:foo\tSM:bar'. [null]")

        # [Options for 'bwa sampe':]
        self.add_option('sampe-a', int, optional=True,
                        description="Maximum insert size for a read pair to be"
                        " considered being mapped properly. Since 0.4.5, this "
                        "option is only used when there are not enough good "
                        "alignment to infer the distribution of insert sizes. "
                        "[500]")
        self.add_option(
            'sampe-o',
            int,
            optional=True,
            description="Maximum occurrences of a read for "
            "pairing. A read with more occurrneces will be treated "
            "as a single-end read. Reducing this parameter helps "
            "faster pairing. [100000]")
        self.add_option('sampe-P', bool, optional=True,
                        description="Load the entire FM-index into memory to "
                        "reduce disk operations (base-space reads only). With "
                        "this option, at least 1.25N bytes of memory are "
                        "required, where N is the length of the genome.")
        self.add_option('sampe-n', int, optional=True,
                        description="Maximum number of alignments to output "
                        "in the XA tag for reads paired properly. If a read "
                        "has more than INT hits, the XA tag will not be "
                        "written. [3]")
        self.add_option(
            'sampe-N',
            int,
            optional=True,
            description="Maximum number of alignments to output "
            "in the XA tag for disconcordant read pairs (excluding "
            "singletons). If a read has more than INT hits, the XA "
            "tag will not be written. [10]")
        self.add_option('sampe-r', str, optional=True,
                        description="Specify the read group in a format like "
                        "'@RG\tID:foo\tSM:bar'. [null]")

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional=True, default="2M")

        # [Options for 'pigz':]
        self.add_option('pigz-blocksize', str, optional=True,
                        default="2048")

    def runs(self, run_ids_connections_files):

        # Check if index is valid
        if not os.path.exists(self.get_option('index') + '.bwt'):
            raise StepError(self, "Could not find index: %s.*" %
                            self.get_option('index'))
        # Compile the list of options
        options_bwa_aln = [
            'aln-n',
            'aln-o',
            'aln-e',
            'aln-d',
            'aln-i',
            'aln-l',
            'aln-k',
            'aln-t',
            'aln-M',
            'aln-E',
            'aln-R',
            'aln-c',
            'aln-N',
            'aln-q',
            'aln-I',
            'aln-B',
            'aln-b',
            'aln-0',
            'aln-1',
            'aln-2']
        options_bwa_samse = ['samse-n', 'samse-r']
        options_bwa_sampe = ['sampe-a', 'sampe-o', 'sampe-P', 'sampe-n',
                             'sampe-N', 'sampe-r']

        set_bwa_aln_options = [o for o in options_bwa_aln if
                               self.is_option_set_in_config(o)]
        set_bwa_samse_options = [o for o in options_bwa_samse if
                                 self.is_option_set_in_config(o)]
        set_bwa_sampe_options = [o for o in options_bwa_sampe if
                                 self.is_option_set_in_config(o)]

        def make_option_list(set_options, prefix=""):
            option_list = list()
            for option in set_options:
                option_wo_prefix = re.sub(r'^%s' % prefix, '', option)
                if isinstance(self.get_option(option), bool):
                    if self.get_option(option):
                        option_list.append('-%s' % option_wo_prefix)
                else:
                    option_list.append('-%s' % option_wo_prefix)
                    option_list.append(str(self.get_option(option)))
            return option_list

        option_list_bwa_aln = make_option_list(set_bwa_aln_options,
                                               prefix="aln-")

        # aln-t option can overwrite default # of cores for bwa aln
        # and the cores variable
        if 'aln-t' not in option_list_bwa_aln:
            option_list_bwa_aln.append('-t')
            option_list_bwa_aln.append(str(self.get_cores()))
        else:
            self.set_cores(self.get_option('aln-t'))

        option_list_bwa_samse = make_option_list(set_bwa_samse_options,
                                                 prefix="samse-")
        option_list_bwa_sampe = make_option_list(set_bwa_sampe_options,
                                                 prefix="sampe-")

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                input_paths = [y for x in [fr_input, sr_input]
                               for y in x if y is not None]

                # Do we have paired end data and is it exactly one ?
                is_paired_end = False if sr_input == [None] else True

                # Fail if we don't have exactly one first read file or
                # an empty connection
                if len(fr_input) != 1 or fr_input == [None]:
                    raise StepError(
                        self, "Expected single input file for first read.")
                # Fail if we don't have exactly one second read file in case of
                # paired end reads
                if is_paired_end and len(sr_input) != 1:
                    raise StepError(
                        self, "Expected single input file for second read.")
                input_paths = fr_input  # single element list
                if is_paired_end:
                    input_paths.extend(sr_input)

                # Check file endings for proper type
                for input_path in input_paths:
                    if len([_ for _ in ['fastq', 'fq', 'fq.gz', 'fastq.gz']
                            if input_path.endswith(_)]) != 1:
                        raise StepError(
                            self, "%s possess unknown suffix. "
                            "(None of: fastq, fq, fq.gz, fastq.gz)")
                # BWA can handle only single files for first and second read
                # IMPORTANT: BWA handles gzipped as well as not gzipped files

                def prepare_input(input_path, exec_group):
                    # 1. Create temporary fifo
                    temp_fifo = run.add_temporary_file(
                        'in-fifo-%s' %
                        os.path.basename(input_path))
                    mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                    exec_group.add_command(mkfifo)
                    # 2. Output input to FIFO
                    dd = [self.get_tool('dd'),
                          'bs=%s' % self.get_option('dd-blocksize'),
                          'if=%s' % input_path,
                          'of=%s' % temp_fifo]
                    exec_group.add_command(dd)
                    return (exec_group, temp_fifo)

                # We need to execute 'bwa aln' first, because the output is the
                # input for 'bwa samse/sampe'
                def execute_bwa_aln(input_path):
                    with run.new_exec_group() as exec_group:
                        exec_group, temp_fifo = prepare_input(
                            input_path, exec_group)
                        # 3. Map reads using bwa aln
                        with exec_group.add_pipeline() as bwa_aln_pipe:
                            # 3.1 Assemble bwa aln command
                            bwa_aln = [
                                self.get_tool('bwa'),
                                'aln'
                            ]
                            bwa_aln.extend(option_list_bwa_aln)
                            bwa_aln.append(self.get_option('index'))
                            bwa_aln.append(temp_fifo)
                            # 3.1.1 Add 'bwa aln' to pipeline
                            bwa_aln_pipe.add_command(bwa_aln)

                            # 3.2 Write bwa aln output to temporary file
                            temp_file = run.add_temporary_file(
                                suffix='%s-bwa-aln.sai' % run_id
                            )
                            dd = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize'),
                                'of=%s' % temp_file
                            ]
                            bwa_aln_pipe.add_command(dd)

                        return (temp_file)

                temp_fr_sai, temp_sr_sai = (str(), str())
                temp_fr_sai = execute_bwa_aln(fr_input[0])
                # And if we handle paired end data
                if is_paired_end:
                    temp_sr_sai = execute_bwa_aln(sr_input[0])

                # Convert the created SAI files to SAM
                with run.new_exec_group() as exec_group:
                    fr_sai_fifo, sr_sai_fifo = (str(), str())
                    temp_fr_fifo, temp_sr_fifo = (str(), str())
                    # 1. Prepare input for first read (fr)
                    # 1.1 The SAI file
                    exec_group, fr_sai_fifo = prepare_input(
                        temp_fr_sai, exec_group)
                    # 1.2 The FASTQ file
                    exec_group, temp_fr_fifo = prepare_input(
                        fr_input[0], exec_group)
                    # 2. Prepare input for second read (sr), if need be
                    if is_paired_end:
                        # 2.1 The SAI file
                        exec_group, sr_sai_fifo = prepare_input(
                            temp_sr_sai, exec_group)
                        # 2.2 The FASTQ file
                        exec_group, temp_sr_fifo = prepare_input(
                            sr_input[0], exec_group)
                        with exec_group.add_pipeline() as bwa_sampe_pipe:
                            # 1. Assemble 'bwa sampe' command
                            bwa_sampe = [
                                self.get_tool('bwa'),
                                'sampe'
                            ]
                            bwa_sampe.extend(option_list_bwa_sampe)
                            bwa_sampe.append(self.get_option('index'))
                            bwa_sampe.append(fr_sai_fifo)
                            bwa_sampe.append(sr_sai_fifo)
                            bwa_sampe.append(temp_fr_fifo)
                            bwa_sampe.append(temp_sr_fifo)
                            # 1.1 Add 'bwa sampe' to pipeline
                            bwa_sampe_pipe.add_command(bwa_sampe)

                            # 2. Compress 'bwa sampe' output
                            pigz = [self.get_tool('pigz'),
                                    '--processes',
                                    str(self.get_cores()),
                                    '--blocksize',
                                    self.get_option('pigz-blocksize'),
                                    '--stdout']
                            bwa_sampe_pipe.add_command(pigz)
                            # 3. Write 'bwa sampe' output to file
                            dd = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize'),
                                'of=%s' %
                                run.add_output_file(
                                    'alignments',
                                    '%s-bwa-aln.sam.gz' % run_id,
                                    input_paths
                                )
                            ]
                            bwa_sampe_pipe.add_command(dd)
                    else:
                        with exec_group.add_pipeline() as bwa_samse_pipe:
                            # 1. Assemble 'bwa samse' command
                            bwa_samse = [
                                self.get_tool('bwa'),
                                'samse'
                            ]
                            bwa_samse.extend(option_list_bwa_samse)
                            bwa_samse.append(self.get_option('index'))
                            bwa_samse.append(fr_sai_fifo)
                            bwa_samse.append(temp_fr_fifo)
                            # 1.1 Add 'bwa samse' to pipeline
                            bwa_samse_pipe.add_command(bwa_samse)

                            # 2. Compress 'bwa samse' output
                            pigz = [self.get_tool('pigz'),
                                    '--processes',
                                    str(self.get_cores()),
                                    '--blocksize',
                                    self.get_option('pigz-blocksize'),
                                    '--stdout']
                            bwa_samse_pipe.add_command(pigz)
                            # 3. Write 'bwa samse' output to file
                            dd = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize'),
                                'of=%s' %
                                run.add_output_file(
                                    'alignments',
                                    '%s-bwa-aln.sam.gz' % run_id,
                                    input_paths
                                )
                            ]
                            bwa_samse_pipe.add_command(dd)
