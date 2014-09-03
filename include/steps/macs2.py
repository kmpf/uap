import sys
from abstract_step import *
import pipeline
import re
import process_pool
import yaml

class Macs2(AbstractStep):
    
    def __init__(self, pipeline):
        super(Macs2, self).__init__(pipeline)
        
        self.set_cores(4)

        self.add_connection('in/alignments')
        self.add_connection('out/log')
        self.add_connection('out/peaks')
        self.add_connection('out/summits')
        self.add_connection('out/diagnosis')
        self.add_connection('out/model')
        self.add_connection('out/negative-peaks')
 
        self.require_tool('macs2')
        self.require_tool('cat4m')
        self.require_tool('pigz')

#        self.add_option('treatment', str, list) # bekommen wir ja immer
        self.add_option('control', str, list)
        self.add_option('format', str, default='AUTO',
                        choices=['ELAND', 'ELANDMULTI', 'ELANDMULTIPET', 
                                 'ELANDEXPORT', 'BED', 'SAM', 'BAM', 'BAMPE', 
                                 'BOWTIE'])
        self.add_option('genome_size', str, default='2.7e9')
        self.add_option('read-length', int, optional=True)
        self.add_option('qvalue-cutoff', float, optional=True)
        self.add_option('pvalue-cutoff', float, optional=True)
        self.add_option('small-local', str, optional=True)
        self.add_option('large-local', str, optional=True)
        self.add_option('shift', int, optional=True)
        self.add_option('keep-duplicates', bool, optional=True)
        self.add_option('broad', bool, optional=True)
        # use "broad-cutoff" only in conjuction with "broad"
        self.add_option('broad-cutoff', float, optional=True)
        self.add_option('to-large', bool, optional=True)
        self.add_option('down-sample', bool, optional=True)
        self.add_option('store-bedgraph', bool, optional=True)
        self.add_option('call-summits', bool, optional=True)
        self.add_option('verbose', int, default=0, choices=[0, 1, 2, 3],
                        optional=True)

    def declare_runs(self):
        

        validated_control_samples = list()
        control_samples = dict()

        if self.get_option('control'):
            control_samples.update(self.get_option('control'))

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection(
                'in/alignments'):
            with self.declare_run(run_id) as run:

                run.add_output_file('log', '%s-macs2-log.txt' % run_id, 
                                    input_paths)
                run.add_output_file('peaks', '%s-macs2-peaks.bed' % run_id,
                                    input_paths)
                run.add_output_file('peaks', '%s-macs2-peaks.xls' % run_id, 
                                    input_paths)
                run.add_output_file('summits', '%s-macs2-summits.bed' % run_id, 
                                    input_paths)
                run.add_output_file('model', '%s-macs2-model.r' % run_id, 
                                    input_paths)
                run.add_output_file('negative-peaks', '%s-macs2-negative-peaks.xls'
                                    % run_id, input_paths)

                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                run.add_private_info('treatment_files', input_paths)
                
                control_files = list()
                for control_id, treatment_list in control_samples.iteritems():
                    if run_id in treatment_list:
                        paths = self.get_input_files_for_run_id_and_connection(
                            control_id, 'in/alignments')
                        control_files.extend(paths)

                run.add_private_info('control_files', control_files)
               
                print('Treatment: ' + " ".join(input_paths) )
                print('Control(s): ' + " ".join(control_files) )
    
    def execute(self, run_id, run):

        macs_out_directory = self.get_temporary_path('macs2-out')

        for control_file in run.get_private_info('control_files'):
            with process_pool.ProcessPool(self) as pool:
                macs2 = [self.get_tool('macs2'), 'callpeak', '--treatment']
                if not run.has_private_info('treatment_files'):
                    raise StandardError("No treatment files for %s to analyse with macs2" % run_id)
                macs2.extend( [" ".join(run.get_private_info('treatment_files'))] )
                # if we do have control data use it
                if not control_file == None:
                    macs2.extend(['--control', control_file ])

                macs2.extend([
                    '--format', self.get_option('format'),
                    '--name', run_id
                ])
                print(macs2)

                try:
                    os.mkdir(macs_out_directory)
                except OSError:
                    pass
                    os.chdir(macs_out_directory)

                    pool.launch(macs2, stdout_path = run.get_single_output_file_for_annotation('log') )

                    peaks_files = run.get_output_files_for_annotation_and_tags('peaks', ['xls', 'bed'])
        
                    # Rename MACS2 generated output files so the connection

                try:
                    os.rename(os.path.join(macs_out_directory, 
                                           '%s_peaks.bed' % run_id), 
                              peaks_files['bed'])
                except OSError:
                    raise StandardError('No file: %s' % os.path.join( macs_out_directory, '%s_peaks.bed' % run_id))

                try:
                    os.rename(os.path.join(macs_out_directory,
                                           '%s_peaks.xls' % run_id), 
                              peaks_files['xls'])
                except OSError:
                    raise StandardError('No file: %s' % os.path.join(macs_out_directory, '%s_peaks.xls' % run_id))

                try:
                    os.rename(os.path.join(macs_out_directory, '%s_summits.bed' % run_id), run.get_single_output_file_for_annotation('summits'))
                except OSError:
                    raise StandardError('No file: %s' % os.path.join(macs_out_directory, '%s_summits.bed' % run_id))

                try:
                    os.rename(os.path.join(macs_out_directory, '%s_model.r' % run_id), run.get_single_output_file_for_annotation('model'))
                except OSError:
                    file_name = run.get_single_output_file_for_annotation('model')
                    with file(file_name, 'a'):
                        os.utime(file_name, None)

                try:
                    os.rename(os.path.join(macs_out_directory, '%s_negative_peaks.xls' % run_id), run.get_single_output_file_for_annotation('negative-peaks'))
                except OSError:
                    raise StandardError('No file: %s' % os.path.join(macs_out_directory, '%s_negative_peaks.xls' % run_id))
