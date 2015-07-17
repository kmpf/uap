import sys
import yaml

from abstract_step import *
import misc
import process_pool

class Segemehl(AbstractStep):
    '''
    segemehl is a software to map short sequencer reads to reference genomes. 
    Unlike other methods, segemehl is able to detect not only mismatches but 
    also insertions and deletions. Furthermore, segemehl is not limited to a 
    specific read length and is able to mapprimer- or polyadenylation 
    contaminated reads correctly.

    This step creates at first two FIFOs. The first through which segemehl gets 
    its genome data and the second to which it writes unmapped reads::

       mkfifo genome_fifo unmapped_fifo
       cat4m <genome-fasta> -o genome_fifo

    The executed segemehl command is this::

        segemehl -d genome_fifo -i <genome-index-file> -q <read1-fastq> 
                 [-p <read2-fastq>] -u unmapped_fifo -H 1 -t 11 -s -S -D 0
                 -o /dev/stdout |  pigz --blocksize 4096 --processes 2 -c

    The unmapped reads are saved via these commands::

        cat4m unmapped_fifo | pigz --blocksize 4096 --processes 2 -c > 
        <unmapped-fastq>

    '''

    def __init__(self, pipeline):
        super(Segemehl, self).__init__(pipeline)

        self.set_cores(12)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/log')
        
        self.require_tool('cat')
        self.require_tool('pigz')
        self.require_tool('segemehl')

        self.add_option('genome', str)
        self.add_option('index', str)
        self.add_option('paired_end', bool, default='yes')
        
    def runs(self, run_ids_connections_files):
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            run = self.new_run(run_id)
            first_read_file = 
            if len(run_ids_connections_files[run_id]['in/first_read']) == 1:
                first_read_file = run_ids_connections_files[run_id]\
                                  ['in/first_read']
            else:
                raise StandardError("Expected one input file.")

            if len(run_ids_connections_files[run_id]['in/second_read']) == 1 \
               and run_ids_connections_files[run_id]['in/second_read']) == None:
                second_read_files = run_ids_connections_files[run_id]\
                                    ['in/second_read']
            elif run_ids_connections_files[run_id]['in/second_read']) == None:
                second_read_files = None
            else:
                raise StandardError("Expected one input file.")

            with run.new_exec_group() as exec_group:
                fifo_path_genome = pool.get_temporary_fifo(
                    'segemehl-genome-fifo', 'input')
                fifo_path_unmapped = pool.get_temporary_fifo(
                    'segemehl-unmapped-fifo', 'output')
                cat = [self.get_tool('cat'), self.get_option('genome'),
                       '-o', fifo_path_genome]
                cat_command = exec_group.add_command(cat)
                
                with segemehl_exec_group.new_pipeline() as pipeline:
                    segemehl = [
                        self.get_tool('segemehl'),
                        '-d', fifo_path_genome,
                        '-i', self.get_option('index'),
                        '-q', first_read_path
                    ]
                    if second_read_files:
                        segemehl.extend(['-p', second_read_path])

                    segemehl.extend([
                        '-u', fifo_path_unmapped,
                        '-H', '1',
                        '-t', cores,
                        '-s', '-S',
                        '-D', '0',
                    ])
                    segemehl_command = pipeline.add_command(
                        segemehl,
                        stderr_path = run.add_output_file(
                            'log', '%s-segemehl-log.txt' % run_id, input_paths)
                    )

                    pigz = [
                        self.get_tool('pigz'), 
                        '--blocksize', '4096', 
                        '--processes', '2', '-c'
                    ]

                    pigz_command = pipeline.add_command(
                        pigz,
                        stdout_path = run.add_output_file(
                            'alignments', '%s-segemehl-results.sam.gz' 
                            % run_id, input_paths)
                    )


                run.add_output_file('unmapped', '%s-segemehl-unmapped.fastq.gz'
                                    % run_id, input_paths)

run.get_single_output_file_for_annotation(
                            'log')
                    )

                )
                pipeline.append(
                    pigz, 
                    stdout_path = run.get_single_output_file_for_annotation(
                        'alignments')
                )

                    
                    
#            pool.launch([self.get_tool('cat4m'), self.get_option('genome'), 
#                         '-o', fifo_path_genome])
            
#            with pool.Pipeline(pool) as pipeline:
#                segemehl = [
#                    self.get_tool('segemehl'),
#                    '-d', fifo_path_genome,
#                    '-i', self.get_option('index'),
#                    '-q', first_read_path
#                    ]

#                if is_paired_end:
#                    segemehl.extend(['-p', second_read_path])
                cores = str(self.get_cores() - 4)

#                segemehl.extend([
#                    '-u', fifo_path_unmapped,
#                    '-H', '1',
#                    '-t', cores,
#                    '-s', '-S',
#                    '-D', '0',
#                ])
#                
#                pigz = [
#                    self.get_tool('pigz'), 
#                    '--blocksize', '4096', 
#                    '--processes', '2', '-c'
#                ]
                pipeline.append(
                    segemehl, 
                    stderr_path = run.get_single_output_file_for_annotation(
                        'log')
                )
                pipeline.append(
                    pigz, 
                    stdout_path = run.get_single_output_file_for_annotation(
                        'alignments')
                )
                
            with pool.Pipeline(pool) as pipeline:
                pigz = [
                    self.get_tool('pigz'), 
                    '--blocksize', '4096', 
                    '--processes', '2', '-c'
                ]
                
                pipeline.append([self.get_tool('cat4m'), fifo_path_unmapped])
                pipeline.append(
                    pigz, 
                    stdout_path = run.get_single_output_file_for_annotation(
                        'unmapped')
                )



    def declare_runs(self):
        found_files = dict()
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        paired_end_info = dict()

        for read in read_types.keys():
            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/%s' % read):
                if input_paths != [None]:
                    if not run_id in paired_end_info:
                        paired_end_info[run_id] = dict()
                    paired_end_info[run_id] = self.\
                                              find_upstream_info_for_input_paths(
                                                  input_paths, 'paired_end')

                # save information per file in found_files
                if not run_id in found_files:
                    found_files[run_id] = dict()

                if not read in found_files[run_id]:
                    found_files[run_id][read] = list()
                # Check if we get exactly one input file
                if len(input_paths) != 1:
                    raise StandardError("Expected one input file.")
                found_files[run_id][read].extend(input_paths)

        for run_id in found_files.keys():
            with self.declare_run(run_id) as run:
                run.add_private_info('paired_end', paired_end_info[run_id])
                input_paths = list()
                for read in found_files[run_id].keys():
                    run.add_private_info(read, found_files[run_id][read])
                    input_paths.extend(found_files[run_id][read])
                
                run.add_output_file('alignments', '%s-segemehl-results.sam.gz' 
                                    % run_id, input_paths)
                run.add_output_file('unmapped', '%s-segemehl-unmapped.fastq.gz'
                                    % run_id, input_paths)
                run.add_output_file('log', '%s-segemehl-log.txt' % run_id, 
                                    input_paths)

    def execute(self, run_id, run):
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        is_paired_end = run.get_private_info('paired_end')
        first_read_path = run.get_private_info('first_read')[0]
        second_read_path = None
        if is_paired_end:
            second_read_path = run.get_private_info('second_read')[0]
        with process_pool.ProcessPool(self) as pool:
            fifo_path_genome = pool.get_temporary_fifo(
                'segemehl-genome-fifo', 'input')
            fifo_path_unmapped = pool.get_temporary_fifo(
                'segemehl-unmapped-fifo', 'output')
            
            pool.launch([self.get_tool('cat4m'), self.get_option('genome'), 
                         '-o', fifo_path_genome])
            
            with pool.Pipeline(pool) as pipeline:
                segemehl = [
                    self.get_tool('segemehl'),
                    '-d', fifo_path_genome,
                    '-i', self.get_option('index'),
                    '-q', first_read_path
                    ]

                if is_paired_end:
                    segemehl.extend(['-p', second_read_path])
                cores = str(self.get_cores() - 4)

                segemehl.extend([
                    '-u', fifo_path_unmapped,
                    '-H', '1',
                    '-t', cores,
                    '-s', '-S',
                    '-D', '0',
                ])
                
                pigz = [
                    self.get_tool('pigz'), 
                    '--blocksize', '4096', 
                    '--processes', '2', '-c'
                ]
                pipeline.append(
                    segemehl, 
                    stderr_path = run.get_single_output_file_for_annotation(
                        'log')
                )
                pipeline.append(
                    pigz, 
                    stdout_path = run.get_single_output_file_for_annotation(
                        'alignments')
                )
                
            with pool.Pipeline(pool) as pipeline:
                pigz = [
                    self.get_tool('pigz'), 
                    '--blocksize', '4096', 
                    '--processes', '2', '-c'
                ]
                
                pipeline.append([self.get_tool('cat4m'), fifo_path_unmapped])
                pipeline.append(
                    pigz, 
                    stdout_path = run.get_single_output_file_for_annotation(
                        'unmapped')
                )

