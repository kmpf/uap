import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml


class TopHat2(AbstractStep):
    '''
    TopHat is a fast splice junction mapper for RNA-Seq reads.
    It aligns RNA-Seq reads to mammalian-sized genomes using the ultra
    high-throughput short read aligner Bowtie, and then analyzes the mapping
    results to identify splice junctions between exons.

    http://tophat.cbcb.umd.edu/
    
    typical command line::

        tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]


    '''
    
    def __init__(self, pipeline):
        super(TopHat2, self).__init__(pipeline)
        self.set_cores(6)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/insertions')
        self.add_connection('out/deletions')
        self.add_connection('out/junctions')
        self.add_connection('out/misc_logs')
        self.add_connection('out/log_stderr')
        self.add_connection('out/prep_reads')
        self.add_connection('out/align_summary')

        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('mkdir')
        self.require_tool('mkfifo')
        self.require_tool('mv')
        self.require_tool('pigz')
#        self.require_tool('bowtie2')
        self.require_tool('tar')
        self.require_tool('tophat2')

#        self.add_option('genome', str)
        self.add_option('index', str)

        self.add_option('library_type', str, choices = [
            'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'])
        
        # TODO: remove swap_reads option
        self.add_option('swap_reads', bool, default = False)

    def runs(self, run_ids_connections_files):

        # Check if option values are valid
        if not os.path.exists(self.get_option('index') + '.1.bt2'):
            raise StandardError("Could not find index file: %s.*" %
                                self.get_option('index') )


        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

#                if len(fr_input) != 1:
#                    raise StandardError("Expected single input file. "
#                                        "Found files %s for run: %s" %
#                                        (fr_input, run_id) )

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
                    # 1. Create temporary fifos and dd files for first read
                    for input_path in fr_input:
                        temp_fifo = run.add_temporary_file(
                            'in-fifo-%s' %
                            os.path.basename(input_path) )
                        mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                        fr_temp_fifos.append(temp_fifo)
                        exec_group.add_command(mkfifo)
                        dd = [self.get_tool('dd'),
                                 'bs=4M',
                                 'if=%s' % input_path,
                                 'of=%s' % temp_fifo]
                        exec_group.add_command(dd)
                    # And if we handle paired end data 
                    if is_paired_end:
                        # 2. Create temporary fifos and dd files for second read
                        for input_path in sr_input:
                            temp_fifo = run.add_temporary_file(
                                'in-fifo-%s' %
                                os.path.basename(input_path) )
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            sr_temp_fifos.append(temp_fifo)
                            exec_group.add_command(mkfifo)
                            dd = [self.get_tool('dd'),
                                  'bs=4M',
                                  'if=%s' % input_path,
                                  'of=%s' % temp_fifo]
                            exec_group.add_command(dd)

                    # 2. Create temporary directory for tophat2 output
                    temp_out_dir = run.add_temporary_directory(
                        "tophat-%s" % run_id)
                    mkdir = [self.get_tool('mkdir'), temp_out_dir]
                    exec_group.add_command(mkdir)

                    # 3. Map reads using tophat2
                    tophat2 = [
                        self.get_tool('tophat2'),
                        '--library-type', self.get_option('library_type'),
                        '--output-dir', temp_out_dir,
                        '-p', str(self.get_cores()),
                        self.get_option('index'), 
                        ','.join(fr_temp_fifos)
                    ]
                    
                    if is_paired_end:
                        tophat2.append(','.join(sr_temp_fifos))

                    exec_group.add_command(
                        tophat2,
                        stderr_path = run.add_output_file(
                            'log_stderr', 
                            '%s-tophat2-log_stderr.txt' % run_id, input_paths)
                    )


                    # Move files created by tophat2 to their final location
                    tophat2_generic_files = [
                        'accepted_hits.bam', 'unmapped.bam', 'insertions.bed',
                        'deletions.bed', 'junctions.bed', 'prep_reads.info',
                        'align_summary.txt'
                    ]

                    # Define output files
                    tophat2_files = {
                        'accepted_hits.bam' : run.add_output_file(
                            'alignments', 
                            '%s-tophat2-accepted.bam' % run_id,
                            input_paths),
                        'unmapped.bam' : run.add_output_file(
                            'unmapped', 
                            '%s-tophat2-unmapped.bam' % run_id,
                            input_paths),
                        'insertions.bed' : run.add_output_file(
                            'insertions', 
                            '%s-tophat2-insertions.bed' % run_id,
                            input_paths),
                        'deletions.bed' : run.add_output_file(
                            'deletions', 
                            '%s-tophat2-deletions.bed' % run_id,
                            input_paths),
                        'junctions.bed' : run.add_output_file(
                            'junctions', 
                            '%s-tophat2-junctions.bed' % run_id,
                            input_paths),
                        'prep_reads.info' : run.add_output_file(
                            'prep_reads', 
                            '%s-tophat2-prep_reads.info' % run_id,
                            input_paths),
                        'align_summary.txt' : run.add_output_file(
                            'align_summary', 
                            '%s-tophat2-align_summary.txt' % run_id,
                            input_paths)
                    }

                # Move files from tophat2 temporary output directory to final
                # destination
                with run.new_exec_group() as clean_up_exec_group:
                    for generic_file, final_path in tophat2_files.items():
                        mv = [self.get_tool('mv'),
                              os.path.join(temp_out_dir, generic_file),
                              final_path
                          ]
                        clean_up_exec_group.add_command(mv)


                    tar_logs = [self.get_tool('tar'),
                                run.add_output_file(
                                    'misc_logs',
                                    '%s-tophat2-misc_logs.tar.gz' % run_id,
                                    input_paths)
                                os.path.join(temp_out_dir, 'logs')
                            ]
                    clean_up_exec_group.add_command(tar_logs )
