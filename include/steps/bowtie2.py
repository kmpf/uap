import sys
import yaml
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Bowtie2(AbstractStep):
    '''
    Bowtie2 is an ultrafast and memory-efficient tool for aligning sequencing
    reads to long reference sequences. It is particularly good at aligning reads
    of about 50 up to 100s or 1,000s of characters, and particularly good at
    aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the
    genome with an FM Index to keep its memory footprint small: for the human
    genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports
    gapped, local, and paired-end alignment modes.

    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    
    typical command line::

        bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]
    '''
    
    def __init__(self, pipeline):
        super(Bowtie2, self).__init__(pipeline)
        self.set_cores(6)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/log_stderr')
        self.add_connection('out/metrics')
        self.add_connection('out/unaligned')


        self.require_tool('pigz')
        self.require_tool('bowtie2')

        self.add_option('index', str, optional=False,
                        description="Path to bowtie2 index (not containing file "
                        "suffixes).")
                        
        self.add_option('cores', int, default=6)
        self.add_option('unaligned', bool, default = None, optional=True)

        # Bowtie2 has so many options that I'm avoiding to add them all now,
        # but it might be necessary later on.

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        # Check if option values are valid
        if not os.path.exists(self.get_option('index') + '.1.bt2'):
            logger.error("Could not find index file: %s.*" %
                         self.get_option('index'))
            sys.exit(1)
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read

                fr_input = run_ids_connections_files[run_id]['in/first_read'][0]
                sr_input = run_ids_connections_files[run_id]['in/second_read'][0]

                # Do we have paired end data and is it exactly one ?
                is_paired_end = True
                input_paths = [fr_input]

                if sr_input == None:
                    is_paired_end = False
                else:
                    input_paths.append(sr_input)


                # bowtie2 is run in this exec group
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as bowtie2_pipe:
                        # Assemble bowtie2 command
                        bowtie2 = [
                            self.get_tool('bowtie2'),
                            '-p', str(self.get_option('cores') - 2),
                            '-x', self.get_option('index')
                        ]
                        if is_paired_end:
                            bowtie2.extend([
                                '-1', fr_input,
                                '-2', sr_input])
                        else:
                            bowtie2.extend(['-U', fr_input])
                
                        log_stderr = run.add_output_file(
                                'log_stderr',
                                '%s-bowtie2-log_stderr.txt' % run_id,
                                                               input_paths)


                        metrics = run.add_output_file(
                                'metrics',
                                '%s-bowtie2-metrics.txt' % run_id,
                            input_paths)
                        bowtie2.extend(['--met-file', metrics])

                        if self.is_option_set_in_config('unaligned'):
                            if self.get_option('unaligned'):
                                unaligned = run.add_output_file(
                                    'unaligned',
                                    '%s-bowtie2-unlagined.fastq' % run_id,
                                    input_paths)
                                bowtie2.extend(['--un', unaligned])

                            
                            


                        bowtie2_pipe.add_command(bowtie2, stderr_path=log_stderr)
                        res = run.add_output_file(
                                'alignments',
                                '%s-bowtie2-results.sam.gz' % run_id,
                                input_paths
                            )
                       
                        # Compress bowtie2 output
                        pigz = [self.get_tool('pigz'),
                                '--stdout']
                        bowtie2_pipe.add_command(pigz, stdout_path=res)

