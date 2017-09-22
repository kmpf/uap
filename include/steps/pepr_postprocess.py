import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class PePrPostprocess(AbstractStep):
    '''
    Post processing of peaks called by PePr.
    Attention: Filter criteria are hard coded!
    '''

    def __init__(self, pipeline):
        super(PePrPostprocess, self).__init__(pipeline)
        
        self.set_cores(4)

        # Mapped Reads
        self.add_connection('in/alignments')
        # Peaks for replicate peak calling
        self.add_connection('in/peaks')

        # Post processed peak lists
        self.add_connection('out/passed_peaks')
        self.add_connection('out/failed_peaks')
        
        self.require_tool('pepr-postprocess')

        # Options for PePr
        ## Required options
        self.add_option('chip_vs_input', dict, optional=False,
                        description='A YAML dictionary that contains: '
                        'runID:                        \n'
                        '    rep1: [<List of runIDs>]  \n'
                        '    input1: [<List of runIDs>]')
        self.add_option('file-format', str, optional=False,
                        choices=['bed', 'sam', 'bam', 'sampe', 'bampe'],
                        description='Read file format. Currently support bed, '
                        'sam, bam, sampe (sam paired-end), bampe '
                        '(bam paired-end)')
        ## Optional options
        self.add_option('remove-artefacts', bool, optional=True, default=True)
        self.add_option('narrow-peak-boundary', bool, optional=True,
                        default=False)
                        
    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['remove-artefacts', 'narrow-peak-boundary']

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
            config_to_option = {'rep1': 'chip','inputs1': 'input'}

            # Get the peak files of runs mentioned in chip_vs_input dict
            try:
                in_files['peak'] = run_ids_connections_files[run_id]['in/peaks']
                if in_files['peak'] == None:
                    logger.error("Upstream run %s provides no peaks" % run_id)
                    sys.exit(1)
                elif in_files['peak'] != 1:
                    logger.error("Expected single peak file for run %s got %s"
                                 % (run_id, len(in_files['peak'])))
                    sys.exit(1)
            except KeyError as e:
                logger.error("No run %s or it provides no peaks" % run_id)
                sys.exit(1)

            # Output file name is coupled to input file name
            file_input_peaks = os.path.basename(in_files['peak'][0])
            (file_passed_peaks, file_failed_peaks) = (file_input_peaks,
                                                      file_input_peaks)
            if self.get_option('remove-artefacts') == True:
                file_passed_peaks += '.passed'
                file_failed_peaks += '.failed'
            if self.get_option('narrow-peak') == True:
                file_passed_peaks += '.boundary_refined'
                file_failed_peaks += '.boundary_refined'
                
            # Get the chip/input files of runs mentioned in chip_vs_input dict
            for key, opt in config_to_option.iteritems():
                experiment = chip_vs_input[run_id]
                in_files[opt] = list()
                try:
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
                # Concatenate all paths of peak/chip/input files with comma
                csv = dict()
                csv = { "--%s" % k: ",".join(v) for k, v in in_files.iteritems()}
                csv = [x for k, v in csv for x in [k, v]]
                input_paths = [f for f in in_files.values()]
                with run.new_exec_group() as pepr_exec_group:
                    # 1. Compile the PePr-postprocess command
                    pepr_post = [self.get_tool('pepr-postprocess'),
                                 '--file-format',
                                 self.get_option('file-format')]
                    ## Append file lists
                    pepr_post.extend(csv)
                    ## Add additional options
                    pepr.extend(option_list)
                    ## Add command to exec group
                    pepr_exec_group.add_command(pepr)
                    
                run.add_output_file('passed_peaks',
                                    file_passed_peaks, input_paths)
                run.add_output_file('failed_peaks',
                                    file_failed_peaks, input_paths)
