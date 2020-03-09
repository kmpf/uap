from uaperrors import StepError
import sys
import os
import re
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Cutadapt(AbstractStep):
    '''
    Cutadapt finds and removes adapter sequences, primers, poly-A tails and
    other types of unwanted sequence from your high-throughput sequencing reads.

    https://cutadapt.readthedocs.org/en/stable/

    This step wraps release: cutadpat 1.5

    '''

    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read', optional=True)
        self.add_connection('out/first_read')
        self.add_connection('out/second_read', optional=True)
        self.add_connection('out/log_first_read')
        self.add_connection('out/log_second_read', optional=True)

        # Step was tested for cat (GNU coreutils) release 8.25
        self.require_tool('cat')
        # Step was tested for cutadapt release 1.16
        self.require_tool('cutadapt')
        # Step was tested for dd (coreutils) release 8.25
        self.require_tool('dd')
        self.require_tool('fix_qnames')
        # Step was tested for mkfifo (GNU coreutils) release 8.25
        self.require_tool('mkfifo')
        # Step was tested for pigz release 2.3.1
        self.require_tool('pigz')

        # Options for cutadapt
        # 1. cutadapt Options that influence how the adapters are found:
        self.add_option("adapter-type", str, optional = True, default='-a',
                        choices=['-a', '-g', '-b'],
                        description="The type of the adapter that has been used "
                        "for sequencing. a: adapter ligated to the 3' end; "
                        "b: adapter ligated to the 3' or 5' end (If the adapter "
                        "is found within the read or overlapping the 3' end of "
                        "the read, the behavior is the same as for the -a value. "
                        "If the adapter overlaps the 5' end (beginning of the read), "
                        "the initial portion of the read matching the adapter is "
                        "trimmed, but anything that follows is kept.); g: adapter "
                        "ligated to the 5' end (If the adapter sequence starts with "
                        "the character '^',  the adapter is 'anchored'. An anchored "
                        "adapter must appear in its entirety at the 5' end of the "
                        "read (it is a prefix of the read). A non-anchored adapter "
                        "may appear partially at the 5' end, or it may occur within "
                        "the read. If it is found within a read, the sequence "
                        "preceding the adapter is also trimmed. In all cases, the "
                        "adapter itself is trimmed).")
        self.add_option("adapter-R1", str, optional = True,
                        description="Adapter sequence to be clipped off of the"
                        "first read.")
        self.add_option("adapter-R2", str, optional = True,
                        description="Adapter sequence to be clipped off of the"
                        "second read")
        self.add_option("adapter-file", str, optional = True,
                        description="File containing adapter sequences to be "
                        "clipped off of the reads.")
        self.add_option("use_reverse_complement", bool, default = False,
                        description="The reverse complement of adapter "
                        "sequences 'adapter-R1' and 'adapter-R2' are used for "
                        "adapter clipping.")
        self.add_option("error-rate", float, optional=True,
                        description="Maximum allowed error rate (no. of errors divided "
                        "by the length of the matching region) (default: 0.1)")
        self.add_option("no-indels", bool, optional=True,
                        description="Do not allow indels in the alignments, that is, "
                        "allow only mismatches. This option is currently only "
                        "supported for anchored 5' adapters (adapter-type: \"-g\" and "
                        "adapter-R[1|2]: \"^ADAPTER\") (default: both mismatches and "
                        "indels are allowed)")
        self.add_option("times", int, optional=True,
                        description="Try to remove adapters at most COUNT times. "
                        "Useful when an adapter gets appended multiple times (default: 1).")
        self.add_option("overlap", int, optional=True,
                        description="Minimum overlap length. If the overlap between the "
                        "read and the adapter is shorter than LENGTH, the read is not "
                        "modified. This reduces the no. of bases trimmed purely due to "
                        "short random adapter matches (default: 3).")
        self.add_option("match-read-wildcards", bool, optional=True,
                        description="Allow 'N's in the read as matches to the adapter "
                        "(default: False).")

        # 2. cutadapt Options for filtering of processed reads:
        self.add_option("discard-trimmed", bool, optional=True,
                        description="Discard reads that contain the adapter instead of "
                        "trimming them. Also use -O in order to avoid throwing away too "
                        "many randomly matching reads!")
        self.add_option("discard-untrimmed", bool, optional=True,
                        description="Discard reads that do not contain the adapter.")
        self.add_option("minimum-length", int, optional=True,
                        description="Discard trimmed reads that are shorter than LENGTH. "
                        "Reads that are too short even before adapter removal are also "
                        "discarded. In colorspace, an initial primer is not counted "
                        "(default: 0).")
        self.add_option("maximum-length", int, optional=True,
                        description="Discard trimmed reads that are longer than LENGTH. "
                        "Reads that are too long even before adapter removal are also "
                        "discarded. In colorspace, an initial primer is not counted "
                        "(default: no limit).")
        self.add_option("no-trim", bool, optional=True,
                        description="Match and redirect reads to output/untrimmed-output as "
                        "usual, but don't remove the adapters. (Default: False)")
        self.add_option("mask-adapter", bool, optional=True,
                        description="Mask with 'N' adapter bases instead of trim (default: "
                        "False)")

        # 3. cutadapt Options that influence what gets output to where:
        # options: [--quiet, --output, --paired-output, --info-file, --rest-file,
        #           --wildcard-file, --too-short-output, --too-long-output,
        #           --untrimmed-output, --untrimmed-paired-output]
        # are handled by uap and not accessible to the user

        # 4. cutadapt Additional modifications to the reads:
        self.add_option("cut", int, optional=True,
                        description="Remove bases from the beginning or end of each read. "
                        "If LENGTH is positive, the bases are removed from the beginning "
                        "of each read. If LENGTH is negative, the bases are removed from "
                        "the end of each read.")
        self.add_option("quality-cutoff", int, optional=True,
                        description="Trim low-quality ends from reads before adapter "
                        "removal. The algorithm is the same as the one used by  BWA "
                        "(Subtract CUTOFF from all qualities; compute partial sums from "
                        "all indices to the end of the sequence; cut sequence at the index "
                        "at which the sum is minimal) (default: 0)")
        self.add_option("quality-base", int, optional=True, choices = [33, 64],
                        description="Assume that quality values are encoded as ascii "
                        "(quality + QUALITY_BASE). The default (33) is usually correct, "
                        "except for reads produced by some versions of the Illumina "
                        "pipeline, where this should be set to 64. (Default: 33)")
        self.add_option("prefix", str, optional=True,
                        description="Add this prefix to read names")
        self.add_option("suffix", str, optional=True,
                        description="Add this suffix to read names")
        self.add_option("strip-suffix", str, optional=True,
                        description="Remove this suffix from read names if present. "
                        "Can be given multiple times.")
        self.add_option("colospace", bool, optional=True,
                        description="Colorspace mode: Also trim the color that is "
                        "adjacent to the found adapter.")
        self.add_option("double-encode", bool, optional=True,
                        description="When in color space, double-encode colors (map "
                        "0,1,2,3,4 to A,C,G,T,N).")
        self.add_option("trim-primer", bool, optional=True,
                        description="When in color space, trim primer base and the first "
                        "color (which is the transition to the first nucleotide)")
        self.add_option("strip-f3", bool, optional=True,
                        description="For color space: Strip the _F3 suffix of read names")
        self.add_option("maq", bool, optional=True,
                        description="MAQ-compatible color space output. This enables "
                        "colorspace, double-encode, trim-primer, strip-f3 and suffix:'/1'.")
        self.add_option("bwa", bool, optional=True,
                        description="BWA-compatible color space output. This enables "
                        "colorspace, double-encode, trim-primer, strip-f3 and suffix:'/1'.")
        self.add_option("length-tag", str, optional=True,
                        description="Search for TAG followed by a decimal number in the "
                        "name of the read (description/comment field of the FASTA or "
                        "FASTQ file). Replace the decimal number with the correct length "
                        "of the trimmed read. For example, use --length-tag 'length=' to "
                        "correct fields like 'length=123'.")
        self.add_option("no-zero-cap", bool, optional=True,
                        description="Do not change negative quality values to zero. "
                        "Colorspace quality values of -1 would appear as spaces in the "
                        "output FASTQ file. Since many tools have problems with that, "
                        "negative qualities are converted to zero when trimming colorspace "
                        "data. Use this option to keep negative qualities.")
        self.add_option("zero-cap", bool, optional=True,
                        description="Change negative quality values to zero. This is "
                        "enabled by default when -c/--colorspace is also enabled. Use "
                        "the above option to disable it.")

        # 5. other non-cutadapt options
        self.add_option('fix_qnames', bool, default=False,
                        description="If set to true, only the leftmost string "
                        "without spaces of the QNAME field of the FASTQ data is "
                        "kept. This might be necessary for downstream analysis.")

        self.add_option('dd-blocksize', str, optional = True, default = "2M")
        self.add_option('pigz-blocksize', str, optional = True, default = "2048")

    def runs(self, cc):

        read_types = {'first_read': 'R1'}

        paired_end = cc.connection_exists('in/second_read')
        if not cc.all_runs_have_connection('in/first_read'):
            read_name = '' if paired_end else ' first'
            run_ids = list(cc.get_runs_without_any('in/first_read'))
            if len(run_ids)>5:
                run_ids = run_ids[0:5] + ['...']
            raise StepError(self, '[cutadapt] No%s read passed by runs '
                           '%s.' % (read_name, list(run_ids)))
        if paired_end:
            if not cc.all_runs_have_connection('in/second_read'):
                read_name = ' second'
                run_ids = list(cc.get_runs_without_any('in/second_read'))
                if len(run_ids)>5:
                    run_ids = run_ids[0:5] + ['...']
                raise StepError(self, '[cutadapt] No%s read passed by runs '
                               '%s.' % (read_name, list(run_ids)))
            read_types['second_read'] = 'R2'

        options = ["error-rate", "no-indels", "times",
                   "overlap", "match-read-wildcards", "discard-trimmed",
                   "discard-untrimmed", "minimum-length", "maximum-length",
                   "no-trim", "mask-adapter", "cut", "quality-cutoff", "quality-base",
                   "prefix", "suffix", "strip-suffix", "colospace", "double-encode",
                   "trim-primer", "strip-f3", "maq", "bwa", "length-tag", "no-zero-cap",
                   "zero-cap"]

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


        for run_id in cc.keys():
            run = self.declare_run(run_id)
            for read in read_types:
                connection = 'in/%s' % read
                input_paths = cc[run_id][connection]

                # make sure that adapter-R1/adapter-R2 or adapter-file are
                # correctly set
                # this kind of mutual exclusive option checking is a bit
                # tedious, so we do it here.
                if read == 'second_read':
                    if ( not self.is_option_set_in_config('adapter-R2') and
                         not self.is_option_set_in_config('adapter-file') ):
                        raise StepError(self,
                            "Option 'adapter-R2' or 'adapter-file' "
                            "required because sample %s is paired end!"
                            % run_id)

                if ( self.is_option_set_in_config('adapter-file') and
                     self.is_option_set_in_config('adapter-R1') ):
                    raise StepError(self,
                        "Option 'adapter-R1' and 'adapter-file' "
                        "are both set but are mutually exclusive!")
                if ( not self.is_option_set_in_config('adapter-file') and
                     not self.is_option_set_in_config('adapter-R1') ):
                    raise StepError(self,
                        "Option 'adapter-R1' or 'adapter-file' "
                        "required to call cutadapt for sample %s!"
                        % run_id)
                temp_fifos = list()
                exec_group = run.new_exec_group()
                for input_path in input_paths:
                    # 1. Create temporary fifo for every input file
                    temp_fifo = run.add_temporary_file(
                        "fifo-%s" % os.path.basename(input_path) )
                    temp_fifos.append(temp_fifo)
                    mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                    exec_group.add_command(mkfifo)
                    # 2. Output files to fifo
                    if input_path.endswith('fastq.gz'):
                        with exec_group.add_pipeline() as pigz_pipe:
                            # 2.1 command: Read file in 4MB chunks
                            dd_in = [
                                self.get_tool('dd'),
                                'ibs=%s' %
                                self.get_option('dd-blocksize'),
                                'if=%s' % input_path
                            ]
                            # 2.2 command: Uncompress file to fifo
                            pigz = [self.get_tool('pigz'),
                                    '--decompress',
                                    '--processes', str(self.get_cores()),
                                    '--blocksize', self.get_option('pigz-blocksize'),
                                    '--stdout']
                            # 2.3 command: Write file in 4MB chunks to
                            #              fifo
                            dd_out = [
                                self.get_tool('dd'),
                                'obs=%s' %
                                self.get_option('dd-blocksize'),
                                'of=%s' % temp_fifo
                            ]

                            pigz_pipe.add_command(dd_in)
                            pigz_pipe.add_command(pigz)
                            pigz_pipe.add_command(dd_out)

                    elif input_path.endswith('fastq'):
                        # 2.1 command: Read file in 4MB chunks and
                        #              write to fifo in 4MB chunks
                        dd_in = [
                            self.get_tool('dd'),
                            'bs=%s' % self.get_option('dd-blocksize'),
                            'if=%s' % input_path,
                            'of=%s' % temp_fifo
                        ]
                        exec_group.add_command(dd_in)
                    else:
                        raise StepError(self, "File %s does not end with any "
                                     "expected suffix (fastq.gz or "
                                     "fastq). Please fix that issue.")
                # 3. Read data from fifos
                with exec_group.add_pipeline() as cutadapt_pipe:
                    # 3.1 command: Read from ALL fifos
                    cat = [self.get_tool('cat')]
                    cat.extend(temp_fifos)
                    cutadapt_pipe.add_command(cat)

                    # 3.2 command: Fix qnames if user wants us to
                    if self.get_option('fix_qnames') == True:
                        fix_qnames = [self.get_tool('fix_qnames')]
                        cutadapt_pipe.add_command(fix_qnames)

                    # Let's get the correct adapter sequences or
                    # adapter sequence fasta file
                    adapter = None
                    # Do we have adapter sequences as input?
                    if self.is_option_set_in_config('adapter-%s' \
                                                    % read_types[read]):
                        # Get adapter sequence
                        adapter = self.get_option(
                            'adapter-%s' % read_types[read])

                        # add index to adapter sequence if necessary
                        if '((INDEX))' in adapter:
                            index = self.find_upstream_info_for_input_paths(
                                input_paths,
                                'index-%s' % read_types[read])
                            adapter = adapter.replace('((INDEX))', index)

                        # create reverse complement if necessary
                        if self.get_option('use_reverse_complement'):
                            complements = adapter.maketrans('acgtACGT',
                                                           'tgcaTGCA')
                            adapter = adapter.translate(complements)[::-1]

                        # make sure the adapter is looking good
                        if re.search('^[ACGT]+$', adapter) == None:
                            raise StepError(self, "Unable to come up with a "
                                         "legit-looking adapter: %s"
                                         % adapter)
                    # Or do we have a adapter sequence fasta file?
                    elif self.is_option_set_in_config('adapter-file'):
                        adapter = "file:" + self.get_option(
                            'adapter-file')
                        if not os.path.exists(
                                self.get_option('adapter-file')):
                            raise StepError(self,
                                "File %s containing adapter sequences "
                                "does not exist."
                                % self.get_option('adapter-file'))


                    # 3.3 command: Clip adapters
                    cutadapt = [self.get_tool('cutadapt'),
                                self.get_option('adapter-type'),
                                adapter, '-']
                    cutadapt.extend(option_list)

                    cutadapt_log_file = run.add_output_file(
                            'log_%s' % read,
                            '%s-cutadapt-%s-log.txt'
                            % (run_id, read_types[read]),
                            input_paths)

                    # 3.4 command: Compress output
                    pigz = [self.get_tool('pigz'),
                            '--processes', str(self.get_cores()),
                            '--blocksize', self.get_option('pigz-blocksize'),
                            '--stdout']
                    # 3.5 command: Write to output file in 4MB chunks
                    clipped_fastq_file = run.add_output_file(
                        "%s" % read,
                        "%s_%s.fastq.gz" %
                        (run_id, read_types[read]),
                        input_paths)

                    dd = [
                        self.get_tool('dd'),
                        'obs=%s' % self.get_option('dd-blocksize'),
                        'of=%s' % clipped_fastq_file
                    ]

                    cutadapt_pipe.add_command(cutadapt,
                                              stderr_path =\
                                              cutadapt_log_file)
                    cutadapt_pipe.add_command(pigz)
                    cutadapt_pipe.add_command(dd)
