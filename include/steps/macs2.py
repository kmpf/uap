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
        self.add_connection('out/diagnosis')
        self.add_connection('out/model')
        # here go either narrow or broad peaks
        self.add_connection('out/peaks')
        self.add_connection('out/peaks-xls')
        # just exist for narrow peaks
        self.add_connection('out/summits')
#        self.add_connection('out/negative-peaks')
        # just exist for broad peaks
        self.add_connection('out/gapped-peaks')



        self.require_tool('macs2')
        self.require_tool('cat4m')
        self.require_tool('pigz')

#        self.add_option('treatment', str, list) # bekommen wir ja immer
        self.add_option('control', dict, optional=False)
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
        
        control_samples = self.get_option('control')

        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection(
                'in/alignments'):

            for control_id, treatment_list in control_samples.iteritems():
                if run_id in treatment_list:
                    new_run_id = "%s-%s" % (run_id, control_id)

                    with self.declare_run(new_run_id) as run:

                        # Files which are created by using --broad
                        if self.is_option_set_in_config('broad'):
                            run.add_output_file('peaks', '%s-macs2_broadPeaks.broadPeak'
                                                % new_run_id, input_paths)
                            run.add_output_file('peaks-xls', '%s-macs2-broadPeaks.xls' 
                                                % new_run_id, input_paths)
                            run.add_output_file('gapped-peaks', '%s-macs2_peaks.gappedPeak'
                                                % new_run_id, input_paths)
                        # Files which are created otherwise
                        if not self.is_option_set_in_config('broad'):
                            run.add_output_file('peaks', '%s-macs2-narrowPeaks.narrowPeak'
                                                % new_run_id, input_paths)
                            run.add_output_file('peaks-xls', '%s-macs2-narrowPeaks.xls' 
                                                % new_run_id, input_paths)
                            run.add_output_file('summits', '%s-macs2-summits.bed' 
                                                % new_run_id, input_paths)

                        run.add_output_file('model', '%s-macs2-model.r' 
                                            % new_run_id, input_paths)
                        run.add_output_file('log', '%s-macs2-log.txt'
                                            % new_run_id, input_paths)
#                        run.add_output_file('negative-peaks', 
#                                            '%s-macs2-negative-peaks.xls'
#                                            % new_run_id, input_paths)

                        if not input_paths:
                            raise StandardError("No input files for run %s" 
                                                % (run_id))
                        run.add_private_info('treatment_files', input_paths)
                
                        control_files =  self.get_input_files_for_run_id_and_connection(control_id, 'in/alignments')

                        run.add_private_info('control_files', control_files)
               
    
    def execute(self, run_id, run):
        # Get the name for a temporary directory
        macs_out_directory = self.get_temporary_path('macs2-out')

        with process_pool.ProcessPool(self) as pool:
            # Assemble MACS2 command
            macs2 = [self.get_tool('macs2'), 'callpeak', '--treatment']
            
            # Fail if there is no treatment file
            if not run.has_private_info('treatment_files'):
                raise StandardError("No treatment files for %s to analyse with macs2" % run_id)
            macs2.extend( [" ".join(run.get_private_info('treatment_files'))] )

            # and if there is no control file
            if not run.has_private_info('control_files'):
                raise StandardError("No control files for %s to analyse with macs2" % run_id)
            macs2.extend(['--control', 
                          " ".join(run.get_private_info('control_files')) ])

            if self.is_option_set_in_config('broad'):
                macs2.extend(['--broad'])

            macs2.extend([
                '--format', self.get_option('format'),
                '--name', run_id,
                '--outdir', macs_out_directory
            ])
            
            try:
                os.mkdir(macs_out_directory)
            except OSError:
                pass
                os.chdir(macs_out_directory)
                
            pool.launch(macs2, stdout_path = 
                        run.get_single_output_file_for_annotation('log') )
                
            # Rename MACS2 generated output files so they can be properly moved
            
        # Define connection:path_to_file dict for general files
        output_files = {'model': os.path.join(macs_out_directory, '%s_model.r' % run_id),
                        'peaks-xls': os.path.join(macs_out_directory, '%s_peaks.xls' % run_id)}

        # Extend dict for broad peaks specific files
        if self.is_option_set_in_config('broad'):
            output_files.update(
                {'peaks': os.path.join(macs_out_directory, '%s_peaks.broadPeak' % run_id),
                 'gapped-peaks': os.path.join(macs_out_directory, '%s_peaks.gappedPeak' % run_id)}
             )

        # Extend dict for narrow peaks specific files
        if not self.is_option_set_in_config('broad'):
            output_files.update(
                {'peaks': os.path.join(macs_out_directory, '%s_peaks.narrowPeak' % run_id),
                 'summits': os.path.join(macs_out_directory, '%s_summits.bed' % run_id)}
             )

        for key, temp_file in output_files.iteritems():
            try:
                os.rename( output_files[key], 
                           run.get_single_output_file_for_annotation(key))
            except OSError:
                raise StandardError('No file: %s' % output_files[key])

        try:
            os.rmdir(macs_out_directory)
        except OSError:
            print("(Should be logged): Directory %s could not be removed. Still files in directory?" 
                  % macs_out_directory)
