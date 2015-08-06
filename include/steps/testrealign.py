import sys
from abstract_step import *
import process_pool
import yaml


class TestRealign(AbstractStep):

    def __init__(self, pipeline):
        super(TestRealign, self).__init__(pipeline)
        
        self.set_cores(12)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/splicesites')
        self.add_connection('out/transrealigned')
        self.add_connection('out/log_stderr')

        self.add_option('genome', str)
        self.add_option('maxdist', int , default=100)
        self.require_tool('testrealign')
        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cat')
        
        






        """
        [-Evn] -d <file> [<file> ...] -q <file> [<file> ...] [-t <n>] [-U <file>] [-T <file>] [-o <file>] [-M <n>]
        Heuristic mapping of short sequences
        time  /home/hackermu/Downloads/segemehl_svn/segemehl/testrealign.x -d /data/bioinf/Data/Indices/hg19/hg19.fa -q bla.sam.gz  -t 12 -o temp
        -d, --database <file> [<file> ...]  list of path/filename(s) of database sequence(s)
        -q, --query <file> [<file> ...]     path/filename of alignment file
        -E, --expand                        expand
        -v, --verbose                       verbose
        -n, --norealign                     do not realign
        -t, --threads <n>                   start <n> threads for realigning (default:1)

        -o, --outfile <file>                path/filename of output sam file (default:none)

        
        """

        


    def declare_runs(self):
        
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                run.add_private_info('in-alignment', input_paths[0])
                run.add_output_file('alignments', '%s-segemehl-realigned.bam' % run_id,input_paths)
                run.add_output_file('splicesites', '%s-segemehl-splicesites.bed' % run_id,input_paths)
                run.add_output_file('transrealigned', '%s-segemehl-transrealigned.bed' % run_id,input_paths)
                run.add_output_file("log_stderr", "%s-segemehl-log_stderr.txt" % run_id, input_paths)
                run.new_exec_group()



    def execute(self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            print('huhu')
            fifo_path_genome = pool.get_temporary_fifo('genome-fifo', 'input')

            
            pool.launch([self.get_tool('cat'), self.get_option('genome'), '-o', fifo_path_genome])
            
            with pool.Pipeline(pool) as pipeline:
                in_alignment = run.get_private_info('in-alignment')                 
                samtools_front = [self.get_tool('samtools'), 'view', '-h', in_alignment]
                
                testrealign = [
                    self.get_tool('testrealign'),
                    '-d', fifo_path_genome,
                    '-q', '/dev/stdin',
                    '-U', run.get_single_output_file_for_annotation('splicesites'),
                    '-T', run.get_single_output_file_for_annotation('transrealigned'),
                    '-t', '6'
                    ]

                samtools_out = [self.get_tool('samtools'), 'view', '-Shbo',   run.get_single_output_file_for_annotation('alignments'), '-']
                pipeline.append(samtools_front)
                pipeline.append(testrealign, stderr_path = run.get_single_output_file_for_annotation('log_stderr'))
                pipeline.append(samtools_out)

            print('huhu')
