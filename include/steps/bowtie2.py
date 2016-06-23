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


        self.require_tool('dd')
        self.require_tool('mkfifo')
        self.require_tool('pigz')
        self.require_tool('bowtie2')

        self.add_option('index', str, optional=False,
                        description="Path to bowtie2 index (not containing file "
                        "suffixes).")
                        
        self.add_option('cores', int, default=6)
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
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                
                input_paths = [ y for x in [fr_input, sr_input] \
                               for y in x if y !=None ]                    

                # Do we have paired end data and is it exactly one ?
                is_paired_end = True
                if sr_input == [None]:
                    is_paired_end = False

                # Tophat is run in this exec group
                with run.new_exec_group() as exec_group:
                    # Lists of fifos
                    fr_temp_fifos = list()
                    sr_temp_fifos = list()
                    # 1. 

                    def prepare_input(input_path, exec_group, temp_fifos):
                        # Create temporary fifo
                        temp_fifo = run.add_temporary_file(
                            'in-fifo-%s' %
                            os.path.basename(input_path), suffix='.fastq' )
                        mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                        exec_group.add_command(mkfifo)
                        temp_fifos.append(temp_fifo)
                        # Is input gzipped fasta?
                        is_fastq_gz = False
                        

                        
                        for suff  in ['fq.gz', 'fastq.gz','.gz']:
                            if input_path.endswith(suff):
                                is_fastq_gz = True

                        # If yes we need to decompress it
                        if is_fastq_gz:
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                      'bs=4M',
                                      'if=%s' % input_path]
                                unzip_pipe.add_command(dd_in)
                                # 2.2 command: Uncompress data
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--stdout']
                                unzip_pipe.add_command(pigz)
                                # 2.3 Write file in 4MB chunks to fifo
                                dd_out = [self.get_tool('dd'),
                                          'obs=4M',
                                          'of=%s' % temp_fifo]
                                unzip_pipe.add_command(dd_out)
                        else:
                            dd = [self.get_tool('dd'),
                                  'bs=4M',
                                  'if=%s' % input_path,
                                  'of=%s' % temp_fifo]
                            exec_group.add_command(dd)

                        return (exec_group, temp_fifos)

                    for input_path in fr_input:
                        exec_group, fr_temp_fifos = prepare_input(
                            input_path, exec_group, fr_temp_fifos)
                    # And if we handle paired end data 
                    if is_paired_end:
                        for input_path in sr_input:
                            exec_group, sr_temp_fifos = prepare_input(
                                input_path, exec_group, sr_temp_fifos)

                    # 3. Map reads using bowtie2
                    with exec_group.add_pipeline() as bowtie2_pipe:
                        # Assemble bowtie2 command
                        bowtie2 = [
                            self.get_tool('bowtie2'), '-q',
                            '-p', str(self.get_option('cores') - 3),
                            '-x', self.get_option('index')
                        ]
                        if is_paired_end:
                            bowtie2.extend([
                                '-1', ','.join(fr_temp_fifos),
                                '-2', ','.join(sr_temp_fifos)])
                        else:
                            bowtie2.extend(['-U', ','.join(fr_temp_fifos)])
                
                        log_stderr = run.add_output_file(
                                'log_stderr',
                                '%s-bowtie2-log_stderr.txt' % run_id,
                                input_paths)

                        metrics = run.add_output_file(
                                'metrics',
                                '%s-bowtie2-metrics.txt' % run_id,
                                input_paths)
                        bowtie2.extend(['--met-file', metrics])
                        bowtie2_pipe.add_command(bowtie2, stderr_path=log_stderr)
                        
                        # Compress bowtie2 output
                        pigz = [self.get_tool('pigz'),
                                '--stdout']
                        bowtie2_pipe.add_command(pigz)
                        # Write bowtie2 output to file
                        dd = [
                            self.get_tool('dd'),
                            'obs=4M',
                            'of=%s' %
                            run.add_output_file(
                                'alignments',
                                '%s-bowtie2-results.sam.gz' % run_id,
                                input_paths
                            )
                        ]
                        bowtie2_pipe.add_command(dd)
