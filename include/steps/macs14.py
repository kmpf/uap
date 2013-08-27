import sys
from abstract_step import *
import pipeline
import re
import process_pool
import yaml

class Macs14(AbstractStep):
    
    def __init__(self, pipeline):
        super(Macs14, self).__init__(pipeline)
        
        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/log')
        self.add_connection('out/peaks')
        self.add_connection('out/summits')
        self.add_connection('out/diagnosis')
#        self.add_connection('out/model')
        self.add_connection('out/negative-peaks')
 
        self.require_tool('macs14')
        self.require_tool('cat4m')
        self.require_tool('pigz')

        self.add_option('control', list , default=list())
        self.add_option('format', str, default='AUTO',
                        choices=['ELAND', 'ELANDMULTI', 'ELANDMULTIPET',
                                 'ELANDEXPORT', 'SAM', 'BAM', 'BOWTIE'])
        self.add_option('genome_size', str, default='hs')


    def declare_runs(self):
        
        control_files = list()
        validated_control_samples = list()
        control_samples = list()

        if self.get_option('control'):
            control_samples.extend(self.get_option('control'))

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:

                run.add_output_file('log', '%s-macs14-log.txt' % run_id, input_paths)
                run.add_output_file('peaks', '%s-macs14-peaks.bed' % run_id, input_paths)
                run.add_output_file('peaks', '%s-macs14-peaks.xls' % run_id, input_paths)
                run.add_output_file('summits', '%s-macs14-summits.bed' % run_id, input_paths)
                run.add_output_file('negative-peaks', '%s-macs14-negative-peaks.xls' 
                                    % run_id, input_paths)

                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                run.add_private_info('treatment_files', input_paths)

                if run_id in control_samples:
                    validated_control_samples.append(run_id)
                    control_files.extend(input_paths)

        for run_id in self.get_run_ids():
            run = self.get_run(run_id)
            if run != None:
                if validated_control_samples and control_files:
                    run.add_private_info('control_files', control_files)
                    run.add_public_info('control_samples', validated_control_samples)
               
    
    def execute(self, run_id, run):

        macs_out_directory = self.get_temporary_path('macs14-out')

        with process_pool.ProcessPool(self) as pool:
            macs14 = [self.get_tool('macs14'), '--treatment']
            if not run.has_private_info('treatment_files'):
                raise StandardError("No treatment files for %s to analyse with macs14" % run_id)
            macs14.extend( [" ".join(run.get_private_info('treatment_files'))] )
            # if we do have control data use it
            if run.has_private_info('control_files'):
                macs14.extend(['--control', " ".join(run.get_private_info('control_files'))])

            macs14.extend([
                    '--format', self.get_option('format'),
                    '--name', run_id
                    ])
            print(macs14)
            try:
                os.mkdir(macs_out_directory)
            except OSError:
                pass
            os.chdir(macs_out_directory)

            pool.launch(macs14, stdout_path = run.get_single_output_file_for_annotation('log') )


        peaks_files = run.get_output_files_for_annotation_and_tags('peaks', ['xls', 'bed'])
        try:
            os.rename(os.path.join(macs_out_directory, '%s_peaks.bed' % run_id), 
                      peaks_files['bed'])
        except OSError:
            raise StandardError('No file: %s' % os.path.join(macs_out_directory, 
                                                             '%s_peaks.bed' % run_id))
        try:
            os.rename(os.path.join(macs_out_directory, '%s_peaks.xls' % run_id), 
                      peaks_files['xls'])
        except OSError:
            raise StandardError('No file: %s' % os.path.join(macs_out_directory, 
                                                             '%s_peaks.xls' % run_id))
        try:
            os.rename(os.path.join(macs_out_directory, '%s_summits.bed' % run_id), 
                      run.get_single_output_file_for_annotation('summits'))
        except OSError:
            raise StandardError('No file: %s' % os.path.join(macs_out_directory, 
                                                             '%s_summits.bed' % run_id))
        try:
            os.rename(os.path.join(macs_out_directory, '%s_negative_peaks.xls' % run_id), 
                      run.get_single_output_file_for_annotation('negative-peaks'))
        except OSError:
            raise StandardError('No file: %s' % os.path.join(macs_out_directory, 
                                                             '%s_negative_peaks.xls' % run_id))
