import sys
from abstract_step import *
import pipeline
import re
import unix_pipeline
import yaml

class Cutadapt(AbstractStep):
    '''
    The cutadapt step can be used to clip adapter sequences from RNASeq reads.
    
    Options:
    
    - ``adapter-R1`` and ``adapter-R2`` for paired-end reads
    - ``adapter`` for non-paired-end reads
    
    Required tools:
    
    - cutadapt (https://code.google.com/p/cutadapt/)
    - pigz (http://zlib.net/pigz/)
    - cat4m
    
    Any adapter may contain ``((INDEX))`` which will be replaced with every
    sample's index. The resulting adapter is checked for sanity and a
    StandardError is thrown if the adapter looks non-legit.
    '''
    
    cores = 6
    connections = []
    connections.append('in/reads')
    connections.append('out/reads')
    connections.append('out/log')
    
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)

    def setup_runs(self, complete_input_run_info):
        # make sure tools are available
        self.tool('cat4m')
        self.tool('pigz')
        self.tool('cutadapt')

        output_run_info = {}
        for step_name, step_input_info in complete_input_run_info.items():
            for input_run_id, input_run_info in step_input_info.items():
                for in_path in sorted(input_run_info['output_files']['reads'].keys()):
                    suffix = ''
                    which = None
                    if input_run_info['info']['paired_end'] == True:
                        which = input_run_info['info']['read_number'][os.path.basename(in_path)]
                        if not which in ['R1', 'R2']:
                            raise StandardError("Expected R1 and R2 input files, but got this: " + in_path)
                        suffix = '-' + which

                    output_run_id = input_run_id + suffix

                    if not output_run_id in output_run_info:
                        output_run_info[output_run_id] = {
                            'output_files': {},
                            'info': {}
                        }
                        if input_run_info['info']['paired_end'] == True:
                            output_run_info[output_run_id]['info']['read_number'] = which

                    # find adapter
                    adapter = self.options['adapter' + suffix]

                    # insert correct index if necessary
                    if '((INDEX))' in adapter:
                        index = input_run_info['info']['index']
                        adapter = adapter.replace('((INDEX))', index)

                    # make sure the adapter is looking good
                    if re.search('^[ACGT]+$', adapter) == None:
                        raise StandardError("Unable to come up with a legit-looking adapter: " + adapter)
                    output_run_info[output_run_id]['info']['adapter'] = adapter

                    for t in [('reads', input_run_id + '-cutadapt' + suffix + '.fastq.gz'),
                            ('log', input_run_id + '-cutadapt' + suffix + '-log.txt')]:
                        pathkey = t[0]
                        path = t[1]
                        if not pathkey in output_run_info[output_run_id]['output_files']:
                            output_run_info[output_run_id]['output_files'][pathkey] = {}
                        if not path in output_run_info[output_run_id]['output_files'][pathkey]:
                            output_run_info[output_run_id]['output_files'][pathkey][path] = []
                        output_run_info[output_run_id]['output_files'][pathkey][path].append(in_path)

        return output_run_info

    def execute(self, run_id, run_info):
        # basic sanity check
        if len(run_info['output_files']['reads']) != 1:
            raise StandardError("Expected a single output file.")

        # set up processes
        cat4m = [self.tool('cat4m')]
        cat4m.extend(*sorted(run_info['output_files']['reads'].values()))

        pigz1 = [self.tool('pigz'), '--processes', '1', '--decompress', '--stdout']

        cutadapt = [self.tool('cutadapt'), '-a', run_info['info']['adapter'], '-']

        pigz2 = [self.tool('pigz'), '--blocksize', '4096', '--processes', '3', '--stdout']

        # create the pipeline and run it
        p = unix_pipeline.create_pipeline()
        p.append(cat4m)
        p.append(pigz1)
        p.append(cutadapt, stderr = open(run_info['output_files']['log'].keys()[0], 'w'))
        p.append(pigz2, stdout = open(run_info['output_files']['reads'].keys()[0], 'w'))

        unix_pipeline.wait()
