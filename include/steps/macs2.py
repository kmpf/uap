from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Macs2(AbstractStep):
    '''
    Model-based Analysis of ChIP-Seq (MACS) is a algorithm, for the identifcation
    of transcript factor binding sites. MACS captures the influence of genome
    complexity to evaluate the significance of enriched ChIP regions, and MACS
    improves the spatial resolution of binding sites through combining the
    information of both sequencing tag position and orientation. MACS can be
    easily used for ChIP-Seq data alone, or with control sample data to increase
    the specificity.

    https://github.com/taoliu/MACS

    typical command line for single-end data::

        macs2 callpeak --treatment <aligned-reads> [--control <aligned-reads>]
                       --name <run-id> --gsize 2.7e9
    '''

    def __init__(self, pipeline):
        super(Macs2, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/log')
        self.add_connection('out/diagnosis')
        self.add_connection('out/model')

        # Narrow peak information
        self.add_connection('out/narrowpeaks')
        self.add_connection('out/narrowpeaks-xls')
        self.add_connection('out/summits')
        # Broad peak information
        self.add_connection('out/broadpeaks')
        self.add_connection('out/broadpeaks-xls')
        self.add_connection('out/gappedpeaks')

        # Step was tested for macs2 release 2.1.1.20160309
        self.require_tool('macs2')
        # Step was tested for mkdir (GNU coreutils) release 8.25
        self.require_tool('mkdir')
        # Step was tested for mv (GNU coreutils) release 8.25
        self.require_tool('mv')

        # Options for MACS2 callpeak subcommand
        ## Input file arguments:
        self.add_option('control', dict, optional=False,
                        description = "Defines the controls and correspondent "
                        "treatments in a YAML hash. Hash keys are the run IDs "
                        "of the control datasets and hash values are the run "
                        "IDs of the treatment datasets.")
        self.add_option('format', str, default='AUTO',
                        choices=['AUTO', 'ELAND', 'ELANDMULTI', 'ELANDMULTIPET',
                                 'ELANDEXPORT', 'BED', 'BEDPE', 'SAM', 'BAM',
                                 'BAMPE', 'BOWTIE'],
                        description = "Format of tag file, can be 'ELAND', "
                        "'BED', 'ELANDMULTI', 'ELANDEXPORT', 'ELANDMULTIPET' "
                        "(for pair-end tags), 'SAM', 'BAM', 'BOWTIE', 'BAMPE' "
                        "or 'BEDPE'. Default is 'AUTO' which will allow MACS "
                        "to decide the format automatically. 'AUTO' is also "
                        "useful when you combine different formats of files. "
                        "Note that MACS can't detect 'BAMPE' or 'BEDPE' format "
                        "with 'AUTO', and you have to implicitly specify the "
                        "format for 'BAMPE' and 'BEDPE'. For more information "
                        "about the formats see https://github.com/taoliu/MACS/")
        self.add_option('gsize', str, default='2.7e9',
                        description = "PLEASE assign this parameter to fit "
                        "your needs! It's the mappable genome size or effective "
                        "genome size which is defined as the genome size which "
                        "can be sequenced. Because of the repetitive features "
                        "on the chromsomes, the actual mappable genome size "
                        "will be smaller than the original size, about 90% or "
                        "70% of the genome size. The default hs -- 2.7e9 is "
                        "recommended for UCSC human hg18 assembly. Here are "
                        "all precompiled parameters for effective genome size: "
                        "hs:2.7e9; mm:1.87e9; ce:9e7; dm:1.2e8")
        self.add_option('tsize', int, optional = True,
                        description = "The size of sequencing tags. If you "
                        "don't specify it, MACS will try to use the first 10 "
                        "sequences from your input treatment file to determine "
                        "the tag size. Specifying it will override the "
                        "automatically determined tag size.")
        self.add_option('bw', int, optional=True,
                        description = "The band width which is used to scan "
                        "the genome ONLY for model building. You can set this "
                        "parameter as the sonication fragment size expected "
                        "from wet experiment. The previous side effect on the "
                        "peak detection process has been removed. So this "
                        "parameter only affects the model building.")
        self.add_option('qvalue', float, optional=True,
                        description = "The qvalue (minimum FDR) cutoff to call "
                        "significant regions. Default is 0.05. For broad marks, "
                        "you can try 0.05 as cutoff. Q-values are calculated "
                        "from p-values using Benjamini-Hochberg procedure.")
        self.add_option('pvalue', float, optional=True,
                        description = "The pvalue cutoff. If 'pvalue' is "
                        "specified, MACS2 will use pvalue instead of qvalue.")
        self.add_option('mfold', str, optional=True,
                        description = "This parameter is used to select the "
                        "regions within MFOLD range of high-confidence "
                        "enrichment ratio against background to build model. "
                        "The regions must be lower than upper limit, and higher "
                        "than the lower limit of fold enrichment. DEFAULT:5,50 "
                        "means using all regions not too low (>5) and not too "
                        "high (<50) to build paired-peaks model. If MACS can "
                        "not find more than 100 regions to build model, it will "
                        "use the --extsize parameter to continue the peak "
                        "detection ONLY if --fix-bimodal is set.")
        self.add_option('nolambda', bool, optional=True,
                        description = "With this flag on, MACS will use the "
                        "background lambda as local lambda. This means MACS "
                        "will not consider the local bias at peak candidate "
                        "regions.")
        self.add_option('slocal', str, optional=True,
                        description = "'slocal' and 'llocal' control which two "
                        "levels of regions will be checked around the peak "
                        "regions to calculate the maximum lambda as local "
                        "lambda. By default, MACS considers 1000bp for small "
                        "local region(--slocal), and 10000bps for large local "
                        "region(--llocal) which captures the bias from a long "
                        "range effect like an open chromatin domain. You can "
                        "tweak these according to your project. Remember that "
                        "if the region is set too small, a sharp spike in the "
                        "input data may kill the significant peak.")
        self.add_option('llocal', str, optional=True,
                        description = "'slocal' and 'llocal' control which two "
                        "levels of regions will be checked around the peak "
                        "regions to calculate the maximum lambda as local "
                        "lambda. By default, MACS considers 1000bp for small "
                        "local region(--slocal), and 10000bps for large local "
                        "region(--llocal) which captures the bias from a long "
                        "range effect like an open chromatin domain. You can "
                        "tweak these according to your project. Remember that "
                        "if the region is set too small, a sharp spike in the "
                        "input data may kill the significant peak.")
        self.add_option('fix-bimodal', bool, optional = True,
                        description = "Whether turn on the auto paired-peak "
                        "model process. If it's set, when MACS failed to build "
                        "paired model, it will use the nomodel settings, the "
                        "'--extsize' parameter to extend each tags. If set, "
                        "MACS will be terminated if paried-peak model is "
                        "failed.")
        self.add_option('nomodel', bool, optional=True,
                        description = "While on, MACS will bypass building the "
                        "shifting model.")
        self.add_option('extsize', int, optional=True,
                        description = "While '--nomodel' is set, MACS uses this "
                        "parameter to extend reads in 5'->3' direction to "
                        "fix-sized fragments. For example, if the size of "
                        "binding region for your transcription factor is 200 "
                        "bp, and you want to bypass the model building by MACS, "
                        "this parameter can be set as 200. This option is only "
                        "valid when --nomodel is set or when MACS fails to "
                        "build model and --fix-bimodal is on.")
        self.add_option('shift', int, optional=True,
                        decsription = "Note, this is NOT the legacy --shiftsize "
                        "option which is replaced by --extsize! You can set an "
                        "arbitrary shift in bp here. Please Use discretion "
                        "while setting it other than default value (0). When "
                        "--nomodel is set, MACS will use this value to move "
                        "cutting ends (5') then apply --extsize from 5' to 3' "
                        "direction to extend them to fragments. When this value "
                        "is negative, ends will be moved toward 3'->5' "
                        "direction, otherwise 5'->3' direction. Recommended to "
                        "keep it as default 0 for ChIP-Seq datasets, or -1 * "
                        "half of EXTSIZE together with --extsize option for "
                        "detecting enriched cutting loci such as certain "
                        "DNAseI-Seq datasets. Note, you can't set values other "
                        "than 0 if format is BAMPE or BEDPE for paired-end "
                        "data. Default is 0. "
                        "Here are some examples for combining --shift and "
                        "--extsize: "
                        "1. To find enriched cutting sites such as some "
                        "DNAse-Seq datasets. In this case, all 5' ends of "
                        "sequenced reads should be extended in both direction "
                        "to smooth the pileup signals. If the wanted smoothing "
                        "window is 200bps, then use '--nomodel --shift -100 "
                        "--extsize 200'. "
                        "2. For certain nucleosome-seq data, we need to pileup "
                        "the centers of nucleosomes using a half-nucleosome "
                        "size for wavelet analysis (e.g. NPS algorithm). Since "
                        "the DNA wrapped on nucleosome is about 147bps, this "
                        "option can be used: '--nomodel --shift 37 --extsize "
                        "73'.")
        self.add_option('keep-dup', int, optional=True,
                        description = "It controls the MACS behavior towards "
                        "duplicate tags at the exact same location -- the same "
                        "coordination and the same strand. The default 'auto' "
                        "option makes MACS calculate the maximum tags at the "
                        "exact same location based on binomal distribution "
                        "using 1e-5 as pvalue cutoff; and the 'all' option "
                        "keeps every tags. If an integer is given, at most this "
                        "number of tags will be kept at the same location. The "
                        "default is to keep one tag at the same location. "
                        "Default: 1")
        self.add_option('broad', bool, default = False, optional=True,
                        description = "When this flag is on, MACS will try to "
                        "composite broad regions in BED12 ( a gene-model-like "
                        "format ) by putting nearby highly enriched regions "
                        "into a broad region with loose cutoff. The broad "
                        "region is controlled by another cutoff through "
                        "--broad-cutoff. The maximum length of broad region "
                        "length is 4 times of d from MACS. DEFAULT: False")
        # use "broad-cutoff" only in conjuction with "broad"
        self.add_option('broad-cutoff', float, optional=True,
                        description = "Cutoff for broad region. This option "
                        "is not available unless --broad is set. If -p is set, "
                        "this is a pvalue cutoff, otherwise, it's a qvalue "
                        "cutoff. DEFAULT: 0.1")
        self.add_option('to-large', bool, optional=True,
                        description = "When set, linearly scale the smaller "
                        "dataset to the same depth as larger dataset, by "
                        "default, the larger dataset will be scaled towards "
                        "the smaller dataset. Beware, to scale up small data "
                        "would cause more false positives.")
        self.add_option('down-sample', bool, optional=True,
                        description = "When set, random sampling method will "
                        "scale down the bigger sample. By default, MACS uses "
                        "linear scaling. This option will make the results "
                        "unstable and irreproducible since each time, random "
                        "reads would be selected, especially the numbers "
                        "(pileup, pvalue, qvalue) would change. Consider to "
                        "use 'randsample' script before MACS2 runs instead.")
        self.add_option('bdg', bool, optional=True,
                        description = "If this flag is on, MACS will store the "
                        "fragment pileup, control lambda, -log10pvalue and "
                        "-log10qvalue scores in bedGraph files. The bedGraph "
                        "files will be stored in current directory named "
                        "NAME+'_treat_pileup.bdg' for treatment data, "
                        "NAME+'_control_lambda.bdg' for local lambda values "
                        "from control, NAME+'_treat_pvalue.bdg' for Poisson "
                        "pvalue scores (in -log10(pvalue) form), and "
                        "NAME+'_treat_qvalue.bdg' for q-value scores from "
                        "Benjamini-Hochberg-Yekutieli procedure "
                        "<http://en.wikipedia.org/wiki/False_discovery_rate#Dependent_tests>")
        self.add_option('call-summits', bool, optional=True,
                        description = "MACS will now reanalyze the shape of "
                        "signal profile (p or q-score depending on cutoff "
                        "setting) to deconvolve subpeaks within each peak "
                        "called from general procedure. It's highly recommended "
                        "to detect adjacent binding events. While used, the "
                        "output subpeaks of a big peak region will have the "
                        "same peak boundaries, and different scores and peak "
                        "summit positions.")
        self.add_option('verbose', int, default=0, choices=[0, 1, 2, 3],
                        optional=True,
                        description = "If you don't want to see any message "
                        "during the running of MACS, set it to 0. But the "
                        "CRITICAL messages will never be hidden. If you want "
                        "to see rich information like how many peaks are "
                        "called for every chromosome, you can set it to 3 or "
                        "larger than 3.")
        # LEGACY options
        self.add_option('buffer-size', int, optional=True,
                        description = "LEGACY option.")
        self.add_option('read-length', int, optional=True,
                        description = "LEGACY option.")

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['format', 'gsize', 'tsize', 'bw', 'qvalue', 'pvalue',
                   'mfold', 'nolambda', 'slocal', 'llocal', 'fix-bimodal',
                   'nomodel', 'extsize', 'shift', 'keep-dup', 'broad',
                   'broad-cutoff', 'to-large', 'down-sample', 'bdg',
                   'call-summits', 'verbose',
                   # LEGACY options
                   'buffer-size', 'read-length']

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

        control_samples = self.get_option('control')
        for control_id, treatment_list in control_samples.items():
            # Check for existence of control files
            control_files = list()
            if control_id != 'None':
                try:
                    control_files = run_ids_connections_files[control_id]\
                                    ['in/alignments']
                    control_id = "-" + control_id
                except KeyError:
                    raise StepError(self, "No control for ID '%s' found." % control_id)
            else:
                control_id = ""

            # Check for existence of treatment files
            for tr in treatment_list:
                treatments = dict()
                try:
                    treatments[tr] = run_ids_connections_files[tr]\
                                     ['in/alignments']
                except KeyError:
                    raise StepError(self, "No treatment for ID '%s' found." % tr)
                # Assemble rund ID
                run_id = "%s%s" % (tr, control_id)

                # Create list of input files
                input_paths = [f for l in [treatments[tr], control_files]\
                               for f in l]

                with self.declare_run(run_id) as run:
                    # Create empty output connections depending on ...
                    result_files = dict()
                    if not self.is_option_set_in_config( 'nomodel' ):
                        result_files["%s_model.r" % run_id] = run.add_output_file(
                            'model',
                            '%s-macs2-model.r' % run_id,
                            input_paths
                        )

                    if not self.get_option('broad'):
                        # ... if we compute narrow peaks ...
                        run.add_empty_output_connection("broadpeaks")
                        run.add_empty_output_connection("broadpeaks-xls")
                        run.add_empty_output_connection("summits")
                        # Result files for narrow peaks
                        narrow_peak = "%s_peaks.narrowPeak" % run_id
                        result_files[narrow_peak] = run.add_output_file(
                            'narrowpeaks',
                            '%s-macs2-narrowPeaks.narrowPeak' % run_id,
                            input_paths
                        )
                        narrow_peak_xls = "%s_peaks.xls" % run_id
                        result_files[narrow_peak_xls] = run.add_output_file(
                            'narrowpeaks-xls',
                            '%s-macs2-narrowPeaks.xls' % run_id,
                            input_paths
                        )
                        summits = "%s_summits.bed" % run_id
                        result_files[summits] = run.add_output_file(
                            'summits',
                            '%s-macs2-summits.bed' % run_id,
                            input_paths
                        )
                    else:
                        # ... or we compute broad peaks.
                        run.add_empty_output_connection("narrowpeaks")
                        run.add_empty_output_connection("narrowpeaks-xls")
                        run.add_empty_output_connection("gappedpeaks")
                        # Files which are created by using --broad
                        broad_peak = "%s_peaks.broadPeak" % run_id
                        result_files[broad_peak] = run.add_output_file(
                            'broadpeaks',
                            '%s-macs2_broadPeaks.broadPeak' % run_id,
                            input_paths
                        )
                        broad_peak_xls = "%s_peaks.xls" % run_id
                        result_files[broad_peak_xls] = run.add_output_file(
                            'broadpeaks-xls',
                            '%s-macs2-broadPeaks.xls' % run_id,
                            input_paths
                        )
                        gapped_peak = "%s_peaks.gappedPeak" % run_id
                        result_files[gapped_peak] = run.add_output_file(
                            'gappedpeaks',
                            '%s-macs2_peaks.gappedPeak' % run_id,
                            input_paths
                        )

                    # Let's compile our commands
                    temp_dir = str
                    with run.new_exec_group() as macs2_exec_group:
                        # 1. Create temporary directory for MACS output
                        temp_dir = run.add_temporary_directory('macs2-out')
                        mkdir = [self.get_tool('mkdir'), temp_dir]
                        macs2_exec_group.add_command(mkdir)
                        # 2. MACS2 command
                        macs2 = [self.get_tool('macs2'), 'callpeak']
                        macs2.append('--treatment')
                        macs2.extend(treatments[tr])
                        ## Append control information
                        if control_files:
                            macs2.append('--control')
                            macs2.extend(control_files)
                        ## Append known info (--name, --directory)
                        macs2.extend([
                            '--name', run_id,
                            '--outdir', temp_dir
                        ])
                        macs2.extend(option_list)
                        macs2_exec_group.add_command(macs2)

                    with run.new_exec_group() as mv_exec_group:
                        for orig, dest_path in result_files.items():
                            # 3. Move file from temp directory to expected
                            #    position
                            orig_path = os.path.join(temp_dir, orig)
                            mv = [self.get_tool('mv'), orig_path, dest_path]
                            mv_exec_group.add_command(mv)
