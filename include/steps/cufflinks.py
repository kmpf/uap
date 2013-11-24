import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml


class CuffLinks(AbstractStep):
    
    def __init__(self, pipeline):
        super(CuffLinks, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/features')
        self.add_connection('out/skipped')
        self.add_connection('out/genes-fpkm')
        self.add_connection('out/isoforms_fpkm')
        self.add_connection('out/log_stderr')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('cufflinks')

        self.add_option('library_type',str)
        self.add_option('use_mask',str, optional=True)

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:

                run.add_output_file('log_stderr', '%s-cufflinks-log.txt'  % run_id, input_paths)
                run.add_output_file('features', '%s-transcripts.gtf'  % run_id, input_paths)
                run.add_output_file('skipped', '%s-skipped.gtf'  % run_id, input_paths)
                run.add_output_file('genes-fpkm', '%s-genes.fpkm_tracking'  % run_id, input_paths)
                run.add_output_file('isoforms_fpkm', '%s-isoforms.fpkm_tracking'  % run_id, input_paths)
                
                run.add_private_info 
                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))             


                if self.is_option_set_in_config('use_mask'):
                    mask = self.get_option('use_mask')
                    if not os.path.exists(mask):
                        raise StandardError('Maskfile not found %s' % mask )
                    else:
                        run.add_private_info('mask_file', mask)

                run.add_private_info('in-bam', input_paths)
                        

    def execute(self, run_id, run):




        cufflinks_out_path = self.get_temporary_path('cufflinks-out')
        
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:

                cores = str(self.get_cores())
                cufflinks = [
                    self.get_tool('cufflinks'),
                    '-o', cufflinks_out_path,
                    '-p', cores,
                    "--library-type=%s" % self.get_option('library_type')
                    ]
            

                if self.is_option_set_in_config('use_mask'):
                    mask = self.get_option('use_mask')
                    cufflinks.extend(['--mask-file', mask])

                cufflinks.extend(run.get_private_info('in-bam'))

            

                print cufflinks


                log_stderr = run.get_single_output_file_for_annotation('log_stderr')
                print log_stderr
          
            
                add_hints = {'writes': [
                        run.get_single_output_file_for_annotation('skipped'),
                        run.get_single_output_file_for_annotation('features'),
                        run.get_single_output_file_for_annotation('genes-fpkm'),
                        run.get_single_output_file_for_annotation('isoforms_fpkm')
                        ]
                            }

                pipeline.append(cufflinks, stderr_path= log_stderr, hints=add_hints)

            
        os.rename(os.path.join(cufflinks_out_path, 'transcripts.gtf'), run.get_single_output_file_for_annotation('features'))
        os.rename(os.path.join(cufflinks_out_path, 'skipped.gtf'), run.get_single_output_file_for_annotation('skipped'))
        os.rename(os.path.join(cufflinks_out_path, 'genes.fpkm_tracking'), run.get_single_output_file_for_annotation('genes-fpkm'))
        os.rename(os.path.join(cufflinks_out_path, 'isoforms.fpkm_tracking'), run.get_single_output_file_for_annotation('isoforms_fpkm'))
        
        try:
            os.rmdir(cufflinks_out_path)
        except OSError:
            pass
