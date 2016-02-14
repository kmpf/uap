import sys
import os
from abstract_step import AbstractStep

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

        self.require_tool('macs2')
        self.require_tool('pigz')
        self.require_tool('mkdir')
        self.require_tool('mv')

        # Options for MACS2 callpeak subcommand
        ## Input file arguments:
        self.add_option('control', dict, optional=False)
        self.add_option('format', str, default='AUTO',
                        choices=['AUTO', 'ELAND', 'ELANDMULTI', 'ELANDMULTIPET',
                                 'ELANDEXPORT', 'BED', 'SAM', 'BAM', 'BAMPE', 
                                 'BOWTIE'])
        self.add_option('gsize', str, default='2.7e9')
        self.add_option('keep-dup', int, optional=True)
        self.add_option('buffer-size', int, optional=True)
        ## Output arguments:
        self.add_option('verbose', int, default=0, choices=[0, 1, 2, 3],
                        optional=True)
        ## Shifting model arguments:
        self.add_option('read-length', int, optional=True)
        self.add_option('shift', int, optional=True)
        ## Peak calling arguments:
        self.add_option('qvalue', float, optional=True)
        self.add_option('pvalue', float, optional=True)
        self.add_option('to-large', bool, optional=True)
        self.add_option('down-sample', bool, optional=True)
        self.add_option('slocal', str, optional=True)
        self.add_option('llocal', str, optional=True)
        self.add_option('broad', bool, default = False, optional=True)
        # use "broad-cutoff" only in conjuction with "broad"
        self.add_option('broad-cutoff', float, optional=True)
        self.add_option('call-summits', bool, optional=True)

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['format', 'gsize', 'keep-dup', 'buffer-size', 'verbose',
                   'read-length', 'shift', 'qvalue', 'pvalue', 'to-large',
                   'down-sample', 'slocal', 'llocal', 'broad', 'broad-cutoff',
                   'call-summits']

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
        for control_id, treatment_list in control_samples.iteritems():
            # Check for existence of control files
            control_files = list()
            if control_id != 'None':
                try:
                    control_files = run_ids_connections_files[control_id]\
                                    ['in/alignments']
                    control_id = "-" + control_id
                except KeyError:
                    raise StandardError("No control for ID '%s' found."
                                        % control_id)
            else:
                control_id = ""

            # Check for existence of treatment files
            for tr in treatment_list:
                treatments = dict()
                try:
                    treatments[tr] = run_ids_connections_files[tr]\
                                     ['in/alignments']
                except KeyError:
                    raise StandardError("No treatment for ID '%s' found." % tr)

                # Assemble rund ID
                run_id = "%s%s" % (tr, control_id)

                # Create list of input files
                input_paths = [f for l in [treatments[tr], control_files]\
                               for f in l]

                with self.declare_run(run_id) as run:
                    # Create empty output connections depending on ...
                    result_files = dict()
                    result_files["%s_model.r" % run_id] = run.add_output_file(
                        'model',
                        '%s-macs2-model.r' % run_id,
                        input_paths
                    )
                    #result_files["%s" % run_id] = run.add_output_file(
                        #'log',
                        #'%s-macs2-log.txt' % run_id,
                        #input_paths
                    #)

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
                        for orig, dest_path in result_files.iteritems():
                            # 3. Move file from temp directory to expected
                            #    position
                            orig_path = os.path.join(temp_dir, orig)
                            mv = [self.get_tool('mv'), orig_path, dest_path]
                            mv_exec_group.add_command(mv)
