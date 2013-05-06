import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml


class TopHat2(AbstractStep):
    
    def __init__(self, pipeline):
        super(TopHat2, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/reads')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/insertions')
        self.add_connection('out/deletions')
        self.add_connection('out/junctions')
        self.add_connection('out/misc_logs')
        self.add_connection('out/log')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('bowtie2')
        self.require_tool('tophat2')
        
    def setup_runs(self, complete_input_run_info, connection_info):
        if not 'library_type' in self.options:
            raise StandardError("Required option is missing: library_type.")
            
        # make sure files are available
        if not os.path.exists(self.options['index'] + '.1.bt2'):
            raise StandardError("Could not find %s file: %s" % (key, self.options[key]))

        if not 'swap_reads' in self.options:
            self.options['swap_reads'] = False
            
        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for run_id, input_run_info in step_input_info.items():

                output_run_info[run_id] = {}
                
                read_files = misc.assign_strings(input_run_info['output_files']['reads'].keys(), ['R1', 'R2'])
                output_run_info[run_id]['info'] = {}
                output_run_info[run_id]['info']['R1-in'] = read_files['R1']
                output_run_info[run_id]['info']['R2-in'] = read_files['R2']
                output_run_info[run_id]['info']['misc-logs'] = run_id + '-misc-logs.yaml'
                output_run_info[run_id]['info']['prep-reads'] = run_id + '-prep_reads.info'

                output_run_info[run_id]['output_files'] = {}
                output_run_info[run_id]['output_files']['alignments']  = {}
                output_run_info[run_id]['output_files']['alignments'][run_id + '-accepted.bam'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['unmapped']  = {}
                output_run_info[run_id]['output_files']['unmapped'][run_id + '-unmapped.bam'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['insertions']  = {}
                output_run_info[run_id]['output_files']['insertions'][run_id + '-insertions.bed'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['deletions']  = {}
                output_run_info[run_id]['output_files']['deletions'][run_id + '-deletions.bed'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['junctions']  = {}
                output_run_info[run_id]['output_files']['junctions'][run_id + '-junctions.bed'] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['misc_logs']  = {}
                output_run_info[run_id]['output_files']['misc_logs'][output_run_info[run_id]['info']['misc-logs']] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['misc_logs'][output_run_info[run_id]['info']['prep-reads']] = input_run_info['output_files']['reads'].keys()
                output_run_info[run_id]['output_files']['log']  = {}
                output_run_info[run_id]['output_files']['log'][run_id + '-tophat2-log.txt'] = input_run_info['output_files']['reads'].keys()

        return output_run_info

    def execute(self, run_id, run_info):
        
        tophat_out_path = self.get_temporary_path('tophat-out')
        
        with process_pool.ProcessPool(self) as pool:
            
            q = run_info['info']['R1-in']
            p = run_info['info']['R2-in']
            if self.options['swap_reads']:
                q = run_info['info']['R2-in']
                p = run_info['info']['R1-in']
                
            tophat2 = [
                self.tool('tophat2'),
                '--library-type', self.options['library_type'],
                '--output-dir', tophat_out_path,
                '-p', '6', self.options['index'], q, p
            ]

            pool.launch(tophat2, stderr_path = run_info['output_files']['log'].keys()[0], hints = {'writes': [run_info['output_files']['alignments'].keys()[0], run_info['output_files']['unmapped'].keys()[0]]})
            
        os.rename(os.path.join(tophat_out_path, 'accepted_hits.bam'), run_info['output_files']['alignments'].keys()[0])
        os.rename(os.path.join(tophat_out_path, 'unmapped.bam'), run_info['output_files']['unmapped'].keys()[0])
        os.rename(os.path.join(tophat_out_path, 'insertions.bed'), run_info['output_files']['insertions'].keys()[0])
        os.rename(os.path.join(tophat_out_path, 'deletions.bed'), run_info['output_files']['deletions'].keys()[0])
        os.rename(os.path.join(tophat_out_path, 'junctions.bed'), run_info['output_files']['junctions'].keys()[0])
        os.rename(os.path.join(tophat_out_path, 'prep_reads.info'), run_info['info']['prep-reads'])
        
        with open(run_info['info']['misc-logs'], 'w') as f:
            logs = {}
            for path in glob.glob(os.path.join(tophat_out_path, 'logs', '*')):
                logs[os.path.basename(path)] = open(path).read()
                os.unlink(path)
            try:
                os.rmdir(os.path.join(tophat_out_path, 'logs'))
            except OSError:
                pass
            f.write(yaml.dump(logs, default_flow_style = False))

        try:
            os.rmdir(tophat_out_path)
        except OSError:
            pass
