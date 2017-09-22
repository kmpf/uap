import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PePr(AbstractStep):
    '''
        
    '''

    def __init__(self, pipeline):
        super(PePr, self).__init__(pipeline)
        
        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/log')
        self.add_connection('out/parameter')

        # Peaks for replicate peak calling
        self.add_connection('out/peaks')
        # Peaks for differential peak calling
        self.add_connection('out/differential_peaks')

        self.require_tool('pepr')
        self.require_tool('tar')
        self.require_tool('mkdir')
        self.require_tool('mv')

        # Options for PePr
        ## Required options
        self.add_option('chip_vs_input', dict, optional=False,
                        description='A YAML dictionary that contains: '
                        'runID:                        \n'
                        '    rep1: [<List of runIDs>]  \n'
                        '    input1: [<List of runIDs>]\n'
                        '[   rep2: [<List of runIDs>]  \n'
                        '    input2: [<List of runIDs>] ]'
                        'rep2 and input2 are optional and will only be used '
                        'for differential peak calling')
        self.add_option('diff', bool, optional=False, description='Tell PePr '
                        'to perform differential binding analysis or not.')
        self.add_option('file-format', str, optional=False, description='Read '
                        'file format. Currently support bed, sam, bam, sampe '
                        '(sam paired-end), bampe (bam paired-end)')

        ## Optional options
        self.add_option('normalization', str, optional = True,
                        choices=['inter-group', 'intra-group', 'scale', 'no'],
                        description='inter-group, intra-group, scale, or no. '
                        'Default is intra-group for peak-calling and '
                        'inter-group for differential binding analysis. PePr '
                        'is using a modified TMM method to normalize for the '
                        'difference in IP efficiencies between samples (see '
                        'the supplementary methods of the paper). It is making '
                        'an implicit assumption that there is substantial '
                        'overlap of peaks in every sample. However, it is '
                        'sometimes not true between groups (for example, '
                        'between TF ChIP-seq and TF knockout). So for '
                        'differential binding analysis, switch to intra-group '
                        'normalization. scale is simply scaling the reads so '
                        'the total library sizes are the same. no '
                        'normalization will not do normalization.')
        self.add_option('peaktype', str, choices=['sharp', 'broad'],
                        optional = True,
                        description='sharp or broad. Default is broad. PePr '
                        'treats broad peaks (like H3K27me3) and sharp peaks '
                        '(like most transcriptions factors) slightly '
                        'different. Specify this option if you know the '
                        'feature of the peaks.')
        self.add_option('shiftsize', str, optional = True,
                        description='Half the fragment size. '
                        'The number of bases to shift forward and reverse '
                        'strand reads toward each other. If not specified by '
                        'user, PePr will empirically estimate this number from '
                        'the data for each ChIP sample.')
        self.add_option('threshold', float, optional = True,
                        description='p-value cutoff. Default:1e-5.')
        self.add_option('windowsize', str, optional = True,
                        description='Sliding window size. '
                        'If not specified by user, PePr will estimate this by '
                        'calculating the average width of potential peaks. '
                        'The lower and upper bound for PePr estimate is 100bp '
                        'and 1000bp. User provided window size is not '
                        'constrained, but we recommend to stay in this range '
                        '(100-1000bp).')

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['file-format','normalization', 'peaktype', 
                   'shiftsize', 'threshold', 'windowsize']

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

        # Get the essential dictionary with information about the relationship
        # between Input and ChIP samples
        chip_vs_input = self.get_option('chip_vs_input')
        # the highest level keys of the dict are the new runID
        for run_id in chip_vs_input.keys():

            in_files = dict()
            config_to_option = dict()
            # Are we going to perform differential peak calling yes or no?
            if not self.get_option('diff'):
                # If not we only use chip1 and input1
                config_to_option = {'rep1': 'chip1','inputs1': 'input1'}
            else:
                # Else we require chip1+input1 and chip2+input2
                config_to_option = {'rep1': 'chip1', 'inputs1': 'input1',
                                    'rep2': 'chip2', 'inputs2': 'input2'}
            # Check the input from the chip_vs_input dict
            for key, opt in config_to_option.iteritems():
                experiment = chip_vs_input[run_id]
                in_files[opt] = list()
                try:
                    # in_run_id: run ID whose in/alignments files
                    #            are used for pepr's --[chip[12]|input[12]]
                    for in_run_id in experiment[key]:
                        in_files[opt].extend(
                            run_ids_connections_files[in_run_id]\
                            ['in/alignments'])
                        if run_ids_connections_files[in_run_id]\
                           ['in/alignments'] == [None]:
                            logger.error("Upstream run %s provides no "
                                         "alignments for run %s"
                                         % (in_run_id, run_id))
                            sys.exit(1)
                except KeyError as e:
                    logger.error("Required key %s missing in 'chip_vs_input' "
                                 "for run %s" % (key, run_id))
                    sys.exit(1)

            # Create a new run named run_id
            with self.declare_run(run_id) as run:
                # Assemble list of all input files
                input_paths = [f for k in in_files for f in in_files[k]]

                # result_files dict:
                ## keys = temporary file names
                ## values = final file names
                result_files = dict()

                # Is differential peak calling happening?
                if self.get_option('diff'):
                    # If yes we do not get any normal peaks
                    run.add_empty_output_connection("peaks")
                    # but we do get two peak lists with differential peaks
                    chip1_file = '%s__PePr_chip1_peaks.bed' % run_id
                    result_files[chip1_file] = run.add_output_file(
                        'differential_peaks', chip1_file, input_paths)
                    chip2_file = '%s__PePr_chip2_peaks.bed' % run_id
                    result_files[chip2_file] = run.add_output_file(
                        'differential_peaks', chip2_file, input_paths)
                else:
                    # If no we do not get any differential_peaks
                    run.add_empty_output_connection("differential_peaks")
                    # but we do get a peak file
                    peaks_file = '%s__PePr_peaks.bed' % run_id
                    result_files[peaks_file] = run.add_output_file(
                        'peaks', peaks_file, input_paths)

                # parameter file used to run PePr with
                parameter_file = '%s__PePr_parameters.txt' % run_id
                result_files[parameter_file] = run.add_output_file(
                    'parameter', parameter_file, input_paths)

                # temp_dir holds temporary directory path
                temp_dir = str()
                with run.new_exec_group() as pepr_exec_group:
                    # 1. Create temporary directory for PePr output
                    temp_dir = os.path.join(
                        run.get_output_directory_du_jour_placeholder(),
                        'pepr-out')
                    mkdir = [self.get_tool('mkdir'), temp_dir]
                    pepr_exec_group.add_command(mkdir)

                    # 2. Compile the PePr command
                    pepr = [self.get_tool('pepr'),
                            '--output-directory', temp_dir,
                            '--file-format',
                            self.get_option('file-format'),
                            '--name', run_id]
                    ## Add '--[chip[12]|input[12]]' and comma separated list of
                    ## alignment files
                    for opt in in_files.keys():
                        pepr.append('--%s' % opt)
                        pepr.append(','.join(in_files[opt]))

                    if self.get_option('diff'):
                        pepr.append('--diff')
                    ## Add additional options
                    pepr.extend(option_list)
                    pepr_exec_group.add_command(pepr)

                with run.new_exec_group() as mv_exec_group:
                    for orig, dest_path in result_files.iteritems():
                        # 3. Move file from temp directory to expected
                        #    position
                        orig_path = os.path.join(temp_dir, orig)
                        mv = [self.get_tool('mv'), orig_path, dest_path]
                        mv_exec_group.add_command(mv)

                with run.new_exec_group() as tar_exec_group:
                    # 
                    log_file = run.add_output_file(
                        'log', '%s__PePr_debug_log.tar.gz' % run_id, input_paths)
                    # We need to compress the temp directory (which should only
                    # contain the log file) and delete all files in there
                    tar = [self.get_tool('tar'),
                            '--create',
                            '--gzip',
                            '--verbose',
                            '--remove-files',
                            '--file=%s' % log_file,
                            temp_dir
                    ]
                    tar_exec_group.add_command(tar)
