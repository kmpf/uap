import sys
from abstract_step import *
import unix_pipeline
import tempfile
import yaml


class Segemehl(AbstractStep):
    def __init__(self, pipeline):
        super(Segemehl, self).__init__(pipeline)
        self.set_cores(12)

    def setup_runs(self, complete_input_run_info):
        # make sure tools are available
        self.tool('segemehl')
        self.tool('pigz')

        # make sure files are available
        for key in ['genome', 'index']:
            if not os.path.exists(self.options[key]):
                raise StandardError("Could not find " + key + " file: " + self.options[key])

        output_run_info = {}
        for run_id, input_run_info in complete_input_run_info.items():
            output_run_info[run_id] = {}
            output_run_info[run_id]['output_files'] = {}
            output_run_info[run_id]['output_files']['alignments']  = {}
            output_run_info[run_id]['output_files']['alignments'][run_id + '-segemehl-results.sam.gz'] = input_run_info['output_files']['reads'].keys()
            output_run_info[run_id]['output_files']['log']  = {}
            output_run_info[run_id]['output_files']['log'][run_id + '-segemehl-log.txt'] = input_run_info['output_files']['reads'].keys()
            read_files = assign_strings(input_run_info['output_files']['reads'].keys(), ['R1', 'R2'])
            output_run_info[run_id]['info'] = {}
            output_run_info[run_id]['info']['R1-in'] = read_files['R1']
            output_run_info[run_id]['info']['R2-in'] = read_files['R2']

        return output_run_info

    def execute(self, run_id, run_info):
        out_name = run_info['output_files']['alignments'].keys()[0]
        if out_name[-3:] != '.gz':
            raise StandardError("Expected .gz in output file name")

        temp_out_name = out_name[:-3]

        _, fifo_path_genome = tempfile.mkstemp('segemehl-genome-fifo')
        os.close(_)
        os.unlink(fifo_path_genome)
        os.mkfifo(fifo_path_genome)

        #_, fifo_path_r1 = tempfile.mkstemp('segemehl-r1-fifo')
        #os.close(_)
        #os.unlink(fifo_path_r1)
        #os.mkfifo(fifo_path_r1)

        #_, fifo_path_r2 = tempfile.mkstemp('segemehl-r2-fifo')
        #os.close(_)
        #os.unlink(fifo_path_r2)
        #os.mkfifo(fifo_path_r2)

        #_, fifo_path_sam = tempfile.mkstemp('segemehl-sam-fifo')
        #os.close(_)
        #os.unlink(fifo_path_sam)
        #os.mkfifo(fifo_path_sam)

        subprocess.Popen(
            [self.tool('cat4m'), self.options['genome'], '-o', fifo_path_genome],
            preexec_fn = os.setsid)

        #subprocess.Popen(
            #[self.tool('cat4m'), fifo_path_sam, '-o', temp_out_name],
            #preexec_fn = os.setsid)
            
        '''
        subprocess.Popen(
            [self.tool('cat4m'), run_info['info']['R1-in'], '-o', fifo_path_r1],
            preexec_fn = os.setsid)

        subprocess.Popen(
            [self.tool('cat4m'), run_info['info']['R2-in'], '-o', fifo_path_r2],
            preexec_fn = os.setsid)
        '''

        segemehl = [
            self.tool('segemehl'),
            '-d', fifo_path_genome,
            '-i', self.options['index'],
            #'-q', fifo_path_r1,
            #'-p', fifo_path_r2,
            '-q', run_info['info']['R1-in'],
            '-p', run_info['info']['R2-in'],
            '-H', '1',
            '-t', '11',
            '-s', '-S',
            '-o', temp_out_name
        ]
        
        up = unix_pipeline.UnixPipeline()
        up.append(segemehl, 
                  stderr = open(run_info['output_files']['log'].keys()[0], 'w'))
        up.run()
        
        pigz = [self.tool('pigz'), '--blocksize', '4096', '--processes', '11', '-c', temp_out_name]
        up = unix_pipeline.UnixPipeline()
        up.append(pigz, stdout = open(out_name, 'w'))
        up.run()

        os.remove(temp_out_name)
        os.unlink(fifo_path_genome)
        #os.unlink(fifo_path_r1)
        #os.unlink(fifo_path_r2)
        #os.unlink(fifo_path_sam)

'''
-------------------------------------
/home/michael/data/genomes/cre4.fasta
-------------------------------------
READS:
  < 32k  27763x

-------------------------------------------------------------------------------------------------------------------
/home/michael/Desktop/rnaseq-pipeline/out/cutadapt-4ad64293/fix_cutadapt-bf21a9e8/Sample_COPD_363-fixed-R1.fastq.gz
-------------------------------------------------------------------------------------------------------------------
READS:
  < 32k   9741x

-------------------------------------------------------------------------------------------------
/home/michael/Desktop/rnaseq-pipeline/out/temp/temp-njw0ja5o/Sample_COPD_363-segemehl-results.sam
-------------------------------------------------------------------------------------------------
READS:
  128k+      1x
  < 32k     60x
WRITES:
  < 32k     57x

-------------------------------------------------------------------------------------------------------------------
/home/michael/Desktop/rnaseq-pipeline/out/cutadapt-4ad64293/fix_cutadapt-bf21a9e8/Sample_COPD_363-fixed-R2.fastq.gz
-------------------------------------------------------------------------------------------------------------------
READS:
  < 32k   9189x
'''