import sys
from abstract_step import *
import os

from logging import getLogger

logger=getLogger('uap_logger')

class CuffLinks(AbstractStep):

    '''
    CuffLinks is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step applies the cufflinks tool which
    assembles transcriptomes from RNA-Seq data and quantifies their
    expression and produces .gtf files with these annotations.
    For details on cufflinks we refer to the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/

    '''
    def __init__(self, pipeline):
        super(CuffLinks, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/features')
        self.add_connection('out/skipped')
        self.add_connection('out/genes-fpkm')
        self.add_connection('out/isoforms_fpkm')
        self.add_connection('out/log_stderr')
        
        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('cufflinks')

        ## options for cufflinks program
        # [General options:]
        # Option --out-put-dir is set automatically
        self.add_option('num-threads', int, optional = True,
                        description = 'Number of threads used during analysis. '
                        'Default: 1')
        self.add_option('seed', int, optional = True,
                        description = 'Value of random number generator seed. '
                        'Default: 0')
        self.add_option('GTF', bool, optional = True,
                        description = 'Quantitate against reference transcript '
                        'annotations. Use with either RABT or ab initio '
                        'assembly is not supported.')
        self.add_option('GTF-guide', bool, optional=True, 
                        description = 'Tells Cufflinks to use the supplied '
                        'reference annotation a GFF file to guide RABT '
                        'assembly. Reference transcripts will be tiled with '
                        'faux-reads to provide additional information in '
                        'assembly. Output will include all reference '
                        'transcripts as well as any novel genes and isoforms '
                        'that are assembled.')
        self.add_option('mask-file', str, optional =True,
                        description = 'Tells Cufflinks to ignore all reads '
                        'that could have come from transcripts in this GTF '
                        'file. We recommend including any annotated rRNA, '
                        'mitochondrial transcripts other abundant transcripts '
                        'you wish to ignore in your analysis in this file. '
                        'Due to variable efficiency of mRNA enrichment methods '
                        'and rRNA depletion kits, masking these transcripts '
                        'often improves the overall robustness of transcript '
                        'abundance estimates.')
        self.add_option('frag-bias-correct', str, optional = True,
                        description = 'Providing Cufflinks with a multifasta '
                        'file via this option instructs it to run our new bias '
                        'detection and correction algorithm which can '
                        'significantly improve accuracy of transcript '
                        'abundance estimates. Default: NULL')
        self.add_option('multi-read-correct', bool, optional = True,
                        description = 'Tells Cufflinks to do an initial '
                        'estimation procedure to more accurately weight reads '
                        'mapping to multiple locations in the genome. '
                        'Default: FALSE')
        self.add_option('library-type', str, optional = False,
                        choices = ['ff-firststrand', 'fr-firststrand',
                                   'ff-secondstrand', 'fr-secondstrand',
                                   'ff-unstranded', 'fr-unstranded',
                                   'transfrags'],
                        description = 'In cases where Cufflinks cannot '
                        'determine the platform and protocol used to generate '
                        'input reads, you can supply this information '
                        'manually, which will allow Cufflinks to infer source '
                        'strand information with certain protocols. The '
                        'available options are listed below. For paired-end '
                        'data, we currently only support protocols where reads '
                        'are point towards each other.'
                        'Library type: fr-unstranded (default); examples: '
                        'Standard Illumina; description: Reads from the '
                        'left-most end of the fragment (in transcript '
                        'coordinates) map to the transcript strand, and the '
                        'right-most end maps to the opposite strand.'
                        'Library type: fr-firststrand; examples: dUTP, NSR, '
                        'NNSR; description: same as fr-unstranded except we '
                        'enforce the rule that the right-most end of the '
                        'fragment (in transcript coordinates) is the first '
                        'sequenced (or only sequenced for single-end reads). '
                        'Equivalently, it is assumed that only the strand '
                        'generated during first strand synthesis is sequenced.'
                        'Library type: fr-secondstrand; examples: Directional '
                        'Illumina (Ligation), Standard SOLiD; same as '
                        'fr-unstranded except we enforce the rule that the '
                        'left-most end of the fragment (in transcript '
                        'coordinates) is the first sequenced (or only sequenced '
                        'for single-end reads). Equivalently, it is assumed '
                        'that only the strand generated during second strand '
                        'synthesis is sequenced. Default: fr-unstranded')
        self.add_option('library-norm-method', str, choices = ['classic-fpkm'],
                        optional = True,
                        description = 'You can control how library sizes (i.e. '
                        'sequencing depths) are normalized in Cufflinks and '
                        'Cuffdiff. Cuffdiff has several methods that require '
                        'multiple libraries in order to work. Library '
                        'normalization methods supported by Cufflinks work on '
                        'one library at a time. Normalization Method supported '
                        'by Cufflinks: classic-fpkm (Library size factor is '
                        'set to 1 - no scaling applied to FPKM values or '
                        'fragment counts. Default: classic-fpkm')

        # [Advanced Abundance Estimation Options:]
        self.add_option('frag-len-mean', int, optional = True,
                        description = 'This is the expected (mean) fragment '
                        'length. The default is 200bp. Note: Cufflinks now '
                        'learns the fragment length mean for each SAM file, so '
                        'using this option is no longer recommended with '
                        'paired-end reads. Default: 200')
        self.add_option('frag-len-std-dev', int, optional = True,
                        description = 'The standard deviation for the '
                        'distribution on fragment lengths. The default is 80bp. '
                        'Note: Cufflinks now learns the fragment length '
                        'standard deviation for each SAM file, so using this '
                        'option is no longer recommended with paired-end reads. '
                        'Default: 80')
        self.add_option('max-mle-iterations', int, optional = True,
                        description = 'Sets the number of iterations allowed '
                        'during maximum likelihood estimation of abundances. '
                        'Default: 5000')
        self.add_option('compatible-hits-norm', bool, optional = True,
                        description = 'With this option, Cufflinks counts only '
                        'those fragments compatible with some reference '
                        'transcript towards the number of mapped hits used in '
                        'the FPKM denominator. This option can be combined with '
                        '-N/--upper-quartile-norm. Default: FALSE')
        self.add_option('total-hits-norm', bool, optional = True,
                        description = 'With this option, Cufflinks counts all '
                        'fragments, including those not compatible with any '
                        'reference transcript, towards the number of mapped '
                        'hits used in the FPKM denominator. This option can be '
                        'combined with -N/–upper-quartile-norm. Default: TRUE')
        self.add_option('num-frag-count-draws', int, optional = True,
                        description = 'Number of fragment generation samples. '
                        'Default: 100')
        self.add_option('num-frag-assign-draws', int, optional = True,
                        description = 'Number of fragment assignment samples '
                        'per generation. Default: 50')
        self.add_option('max-frag-multihits', str, optional = True, 
                        description = 'Maximum number of alignments allowed '
                        'per fragment. Default: unlim')
        self.add_option('no-effective-length-correction', bool, optional = True, 
                        description = 'Cufflinks will not employ its '
                        '"effective" length normalization to transcript FPKM. '
                        'Default: FALSE')
        self.add_option('no-length-correction', bool, optional = True, 
                        description = 'Cufflinks will not normalize fragment '
                        'counts by transcript length at all. Use this option '
                        'when fragment count is independent of the size of the '
                        'features being quantified (e.g. for small RNA '
                        'libraries, where no fragmentation takes place, or 3 '
                        'prime end sequencing, where sampled RNA fragments are '
                        'all essentially the same length). Experimental option, '
                        'use with caution. Default: FALSE')
        self.add_option('upper-quartile-norm', bool, optional = True,
                        description = 'DEPRECATED! Use --library-norm-method '
                        'With this option, Cufflinks normalizes '
                        'by the upper quartile of the number of fragments '
                        'mapping to individual loci instead of the total number '
                        'of sequenced fragments. This can improve robustness of '
                        'differential expression calls for less abundant genes '
                        'and transcripts.')

        # [Advanced Assembly Options:]
        self.add_option('label', str, optional = True, 
                        description = 'Cufflinks will report transfrags in GTF '
                        'format, with a prefix given by this option. Default: '
                        'CUFF')
        self.add_option('min-isoform-fraction', float, optional = True, 
                        description = 'After calculating isoform abundance for '
                        'a gene, Cufflinks filters out transcripts that it '
                        'believes are very low abundance, because isoforms '
                        'expressed at extremely low levels often cannot '
                        'reliably be assembled, and may even be artifacts of '
                        'incompletely spliced precursors of processed '
                        'transcripts. This parameter is also used to filter '
                        'out introns that have far fewer spliced alignments '
                        'supporting them. The default is 0.1, or 10% of the '
                        'most abundant isoform (the major isoform) of the '
                        'gene. Range: 0.0-1.0. Default: 0.10')
        self.add_option('pre-mrna-fraction', float, optional = True, 
                        description = 'Some RNA-Seq protocols produce a '
                        'significant amount of reads that originate from '
                        'incompletely spliced transcripts, and these reads '
                        'can confound the assembly of fully spliced mRNAs. '
                        'Cufflinks uses this parameter to filter out '
                        'alignments that lie within the intronic intervals '
                        'implied by the spliced alignments. The minimum depth '
                        'of coverage in the intronic region covered by the '
                        'alignment is divided by the number of spliced reads, '
                        'and if the result is lower than this parameter value, '
                        'the intronic alignments are ignored. The default is '
                        '15%. Range: 0.0-1.0. Default: 0.15')
        self.add_option('max-intron-length', int, optional = True,
                        description = 'The maximum intron length. Cufflinks '
                        'will not report transcripts with introns longer than '
                        'this, and will ignore SAM alignments with REF_SKIP '
                        'CIGAR operations longer than this. Default: 300000')
        self.add_option('junc-alpha', float, optional = True, 
                        description = 'The alpha value for the binomial test '
                        'used during false positive spliced alignment '
                        'filtration. Default: 0.001')
        self.add_option('small-anchor-fraction', float, optional = True, 
                        description = 'Spliced reads with less than this '
                        'percent of their length on each side of the junction '
                        'are considered suspicious and are candidates for '
                        'filtering prior to assembly. Default: 0.09')
        self.add_option('min-frags-per-transfrag', int, optional = True, 
                        description = 'Assembled transfrags supported by fewer '
                        'than this many aligned RNA-Seq fragments are not '
                        'reported. Default: 10')
        self.add_option('overhang-tolerance', int, optional = True, 
                        description = 'The number of bp allowed to enter the '
                        'intron of a transcript when determining if a read or '
                        'another transcript is mappable to/compatible with it. '
                        'The default is 8 bp based on the default bowtie/'
                        'TopHat parameters. Default: 8')
        self.add_option('max-bundle-length', int, optional = True, 
                        description = 'Maximum genomic length allowed for a '
                        'given bundle. Default: 3500000')
        self.add_option('max-bundle-frags', int, optional = True, 
                        description = 'Sets the maximum number of fragments a '
                        'locus may have before being skipped. Skipped loci '
                        'are listed in skipped.gtf. Default: 500000')
        self.add_option('min-intron-length', int, optional = True, 
                        description = 'Minimum intron size allowed in genome. '
                        'Default: 50')
        self.add_option('trim-3-avgcov-thresh', int, optional = True, 
                        description = 'Minimum average coverage required to '
                        'attempt 3 prime trimming. Default: 10')
        self.add_option('trim-3-dropoff-frac', float, optional = True, 
                        description = 'The fraction of average coverage below '
                        'which to trim the 3 prime end of an assembled '
                        'transcript. Default: 0.1')
        self.add_option('max-multiread-fraction', float, optional = True, 
                        description = 'The fraction a transfrags supporting '
                        'reads that may be multiply mapped to the genome. A '
                        'transcript composed of more than this fraction will '
                        'not be reported by the assembler. Default: 0.75 '
                        '(75% multireads or more is suppressed).')
        self.add_option('overlap-radius', int, optional = True, 
                        description = 'Transfrags that are separated by less '
                        'than this distance (in bp) get merged together, and '
                        'the gap is filled. Default: 50')

        # [Advanced Reference Annotation Guided Assembly Options:]

        self.add_option('no-faux-reads', bool, optional = True, 
                        description = 'This option disables tiling of the '
                        'reference transcripts with faux reads. Use this if '
                        'you only want to use sequencing reads in assembly but '
                        'do not want to output assembled transcripts that lay '
                        'within reference transcripts. All reference '
                        'transcripts in the input annotation will also be '
                        'included in the output. Default: FALSE')
        self.add_option('3-overhang-tolerance', int, optional = True, 
                        description = 'The number of bp allowed to overhang '
                        'the 3 prime end of a reference transcript when '
                        'determining if an assembled transcript should be '
                        'merged with it (i.e., the assembled transcript is not '
                        'novel). Default: 600')
        self.add_option('intron-overhang-tolerance', int, optional = True, 
                        description = 'The number of bp allowed to enter the '
                        'intron of a reference transcript when determining if '
                        'an assembled transcript should be merged with it '
                        '(i.e., the assembled transcript is not novel). '
                        'Default: 50')

        # [Advanced Program Behavior Options:]
        self.add_option('verbose', bool, optional = True, 
                        description = 'Print lots of status updates and other '
                        'diagnostic information. Default: FALSE')
        # Do not allow using 'quiet' option
        self.add_option('no-update-check', bool, optional = True, 
                        description = 'Turns off the automatic routine that '
                        'contacts the Cufflinks server to check for a more '
                        'recent version. Default: FALSE')



    def runs(self, run_ids_connections_files):
        
        # Compile the list of options
        options = [
            # [General options:]
            'num-threads', 'seed', 'GTF', 'GTF-guide', 'mask-file',
            'frag-bias-correct', 'multi-read-correct', 'library-type',
            'library-norm-method',
            # [Advanced Abundance Estimation Options:]
            'frag-len-mean', 'frag-len-std-dev', 'max-mle-iterations',
            'compatible-hits-norm', 'total-hits-norm', 'num-frag-count-draws',
            'num-frag-assign-draws', 'max-frag-multihits',
            'no-effective-length-correction', 'no-length-correction',
            'upper-quartile-norm',
            # [Advanced Assembly Options:]
            'label', 'min-isoform-fraction', 'pre-mrna-fraction',
            'max-intron-length', 'junc-alpha', 'small-anchor-fraction',
            'min-frags-per-transfrag', 'overhang-tolerance',
            'max-bundle-length', 'max-bundle-frags', 'min-intron-length',
            'trim-3-avgcov-thresh', 'trim-3-dropoff-frac',
            'max-multiread-fraction', 'overlap-radius',
            # [Advanced Reference Annotation Guided Assembly Options:]
            'no-faux-reads', '3-overhang-tolerance',
            'intron-overhang-tolerance',
            # [Advanced Program Behavior Options:]
            'verbose', 'no-update-check']

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

        # 'num-threads' option can overwrite default # of cores for cufflinks
        # and the cores variable
        if 'num-threads' not in set_options:
            option_list.append('--num-threads')
            option_list.append(str(self.get_cores()))
        else:
            self.set_cores(self.get_option('num-threads'))


        for run_id in run_ids_connections_files.keys():

             with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/alignments']
                temp_dir = run.add_temporary_directory('cufflinks-out')
                
                # check, if only a single input file is provided
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file, "
                                        "but got this %s" % input_paths)

                cufflinks = [self.get_tool('cufflinks'),'-o', temp_dir, '-q']
                cufflinks.extend(option_list)
                cufflinks.append(input_paths[0])

                # 3. Announce output files
                #    !!!Keep in mind that run.add_output_file() returns the
                #    expected path of the output file

                result_files = {
                    'transcripts.gtf': run.add_output_file(
                        'features',
                        '%s-transcripts.gtf' % run_id,
                        input_paths
                        ),
                    'skipped.gtf': run.add_output_file(
                        'skipped',
                        '%s-skipped.gtf' % run_id,
                        input_paths
                        ),
                    'genes.fpkm_tracking': run.add_output_file(
                        'genes-fpkm',
                        '%s-genes.fpkm_tracking' % run_id,
                        input_paths
                        ),
                    'isoforms.fpkm_tracking': run.add_output_file(
                        'isoforms_fpkm',
                        '%s-isoforms.fpkm_tracking' % run_id,
                        input_paths
                        )
                    }
             
             with run.new_exec_group() as exec_group:
                 # 1. Create temporary directory for cufflinks output
                 mkdir = [self.get_tool('mkdir'), temp_dir]
                 exec_group.add_command(mkdir)

             # 2. Execute cufflinks
             with run.new_exec_group() as exec_group:
                 exec_group.add_command(
                     cufflinks, stderr_path = run.add_output_file(
                         'log_stderr',
                         '%s-cufflinks-log.txt' % run_id,
                         input_paths)
                 )

             # Files get moved to expected position after cufflinks finished
             with run.new_exec_group() as mv_exec_group:
                 for orig, dest_path in result_files.iteritems():
                     # 3. Rename files 
                     orig_path = os.path.join(temp_dir, orig)
                     mv = [self.get_tool('mv'), orig_path, dest_path]
                     mv_exec_group.add_command(mv)
