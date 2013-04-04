import sys
from abstract_step import *
import subprocess

class FixCutadapt(AbstractStep):
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)
        self.set_cores(9)

    def setup_runs(self, input_run_info_complete):
        # make sure tools are available
        self.tool('pigz')

        output_run_info = {}
        for input_run_id, input_run_info in input_run_info_complete.items():
            if not 'read_number' in input_run_info['info']:
                raise StandardError("fix_cutadapt can only be run on paired-end sequenced samples.")
            new_key = input_run_id.replace('-R1', '').replace('-R2', '')
            if not new_key in output_run_info:
                output_run_info[new_key] = { 'output_files': { 'reads': {} }, 'info': {} }
            output_run_info[new_key]['info'][input_run_info['info']['read_number'] + '-in'] = input_run_info['output_files']['reads'].keys()[0]
            for in_path in sorted(input_run_info['output_files']['reads'].keys()):
                for _ in ['R1', 'R2']:
                    k2 = new_key + '-fixed-' + _ + '.fastq.gz'
                    if not k2 in output_run_info[new_key]['output_files']['reads']:
                        output_run_info[new_key]['output_files']['reads'][k2] = []
                    output_run_info[new_key]['output_files']['reads'][k2].append(in_path)
                    output_run_info[new_key]['info'][_ + '-out'] = k2
        return output_run_info

    def execute(self, run_id, run_info):
        p1 = subprocess.Popen([self.tool('pigz'), "--decompress", "--stdout", "--blocksize", "4096", "--processes", "1", run_info['info']['R1-in']], bufsize = -1, stdout = subprocess.PIPE)
        fin1 = p1.stdout
        p2 = subprocess.Popen([self.tool('pigz'), "--decompress", "--stdout", "--blocksize", "4096", "--processes", "1", run_info['info']['R2-in']], bufsize = -1, stdout = subprocess.PIPE)
        fin2 = p2.stdout

        fout1f = open(run_info['info']['R1-out'], 'w')
        fout2f = open(run_info['info']['R2-out'], 'w')

        p3 = subprocess.Popen([self.tool('pigz'), "--blocksize", "4096", "--processes", "3", "--stdout"], bufsize = -1, stdin = subprocess.PIPE, stdout = fout1f)
        fout1 = p3.stdin
        p4 = subprocess.Popen([self.tool('pigz'), "--blocksize", "4096", "--processes", "3", "--stdout"], bufsize = -1, stdin = subprocess.PIPE, stdout = fout2f)
        fout2 = p4.stdin

        rcount = 0
        wcount = 0
        while True:
            lines1 = []
            lines2 = []
            for _ in range(4):
                lines1.append(fin1.readline())
                lines2.append(fin2.readline())
            if (len(lines1[0]) == 0):
                break
            if not (lines1[1] == "\n" or lines2[1] == "\n"):
                for _ in range(4):
                    fout1.write(lines1[_])
                    fout2.write(lines2[_])
                wcount += 1
            rcount += 1

        fout1.flush()
        fout1.close()
        p3.wait()

        fout2.flush()
        fout2.close()
        p4.wait()

        print("Read " + str(rcount) + " entries, wrote " + str(wcount) + " entries.")