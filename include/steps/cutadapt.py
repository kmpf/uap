import sys
from abstract_step import *
import pipes
import yaml

class Cutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)

        # for every output file, remember which read it was (R1 or R2)
        self.options_for_output_file = {}

    def setup_runs(self, input_run_info):
        output_run_info = {}
        for key, input_files in input_run_info.items():
            for fn in input_files:
                if '_R1_' in fn:
                    k = key + '-R1'
                    if not k in output_run_info:
                        output_run_info[k] = {}
                    destination_file_name = key + '-cutadapt-R1.fastq.gz'
                    if not destination_file_name in output_run_info[k]:
                        output_run_info[k][destination_file_name] = []
                    output_run_info[k][destination_file_name].append(fn)
                    self.options_for_output_file[destination_file_name] = {'read': 'R1'}
                elif '_R2_' in fn:
                    k = key + '-R2'
                    if not k in output_run_info:
                        output_run_info[k] = {}
                    destination_file_name = key + '-cutadapt-R2.fastq.gz'
                    if not destination_file_name in output_run_info[k]:
                        output_run_info[k][destination_file_name] = []
                    output_run_info[k][destination_file_name].append(fn)
                    self.options_for_output_file[destination_file_name] = {'read': 'R2'}
                else:
                    raise StandardError("Expected input files with _R1_ or _R2_.")
        return output_run_info

    def execute(self, run_id, run_info):
        print("executing " + self.get_step_id() + " " + run_id)
        #print(yaml.dump(run_info, default_flow_style = False))
        if len(run_info) != 1:
            raise StandardError("Expected a single output file.")
        basename = os.path.basename(run_info.keys()[0])
        options = self.options_for_output_file[basename]
        if options['read'] == 'R1':
            index = ''
            adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" + index + "ATCTCGTATGCCGTCTTCTGCTTG"
            command = []
            command.append(self.pipeline.config['tools']['cutadapt'])
            command.append('-a')
            command.append(adapter)

            t = pipes.Template()
            t.prepend(self.pipeline.config['tools']['pigz']['path'] +
                " -d -c " +
                ' '.join(run_info.values()[0]), ".-")
            t.append(self.pipeline.config['tools']['cutadapt']['path'] +
                " -a " + adapter + " - ", "--")
            t.append(self.pipeline.config['tools']['pigz']['path'] +
                " --blocksize 4096 --processes 3 -c" , "--")
            exit_code = t.copy('/dev/null', run_info.keys()[0])
            if exit_code != 0:
                raise StandardError("Command did not complete successfully.")

            '''
            #command = "python2.7 ../../tools/cutadapt/cutadapt-1.2.1/bin/cutadapt -a " + adapter +
            " \"" + sourcePath + "\" 2> \"" + logPath + "\" | pigz --processes 3 --blocksize 4096
            -c > \"" + destinationPath + ".tmp\""
            '''
        super(Cutadapt, self).execute(run_id, run_info)
