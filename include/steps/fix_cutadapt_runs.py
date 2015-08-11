import sys
from abstract_step import *
import process_pool

class FixCutadapt(AbstractStep):
    '''
    Remove both reads of a paired end read if one of them has been completely
    removed by cutadapt
    '''
    def __init__(self, pipeline):
        super(FixCutadapt, self).__init__(pipeline)
        
        self.set_cores(6)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        
        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('fix_cutadapt')
        self.require_tool('mkfifo')
        self.require_tool('pigz')

    def runs(self, run_ids_connections_files):
        '''

        '''
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                temp_fifos = dict()
                exec_group = run.new_exec_group()
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]
                    temp_fifos["%s_in" % read] = None
                    temp_fifos["%s_out" % read] = None
                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)

                    elif len(input_paths) != 1:
                        raise StandardError("Expected single input file. "
                                            "Found files %s for run: %s" %
                                            (input_paths, run_id) )
                    else:
                        # 1. Create temporary fifos
                        # 1.1 Input fifo
                        temp_fifos["%s_in" % read] = run.add_temporary_file(
                            "in-fifo-%s" %
                            os.path.basename(input_paths[0]) )
                        mkfifo_in = [self.get_tool('mkfifo'),
                                     temp_fifos["%s_in" % read]]
                        exec_group.add_command(mkfifo_in)
                        # 1.2 Output fifo
                        temp_fifos["%s_out" % read] = run.add_temporary_file(
                            "%s-out-fifo" % read)
                        mkfifo_out = [self.get_tool('mkfifo'),
                                      temp_fifos["%s_out" % read]]
                        exec_group.add_command(mkfifo_out)

                        # 2. Output files to fifo
                        if input_paths[0].endswith('fastq.gz'):
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                         'ibs=4M',
                                         'if=%s' % input_paths[0]]
                                # 2.2 command: Uncompress file to fifo
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--stdout']
                                # 2.3 Write file in 4MB chunks to fifo
                                dd_out = [self.get_tool('dd'),
                                          'obs=4M',
                                          'of=%s' % temp_fifos["%s_in" % read] ]

                                unzip_pipe.add_command(dd_in)
                                unzip_pipe.add_command(pigz)
                                unzip_pipe.add_command(dd_out)
                        elif input_paths[0].endswith('fastq'):
                            # 2.1 command: Read file in 4MB chunks and
                            #              write to fifo in 4MB chunks
                            dd_in = [self.get_tool('dd'),
                                     'bs=4M',
                                     'if=%s' % input_paths[0],
                                     'of=%s' % temp_fifos["%s_in" % read] ]
                            exec_group.add_command(dd_in)
                        else:
                            raise StandardError("File %s does not end with "
                                                "any expected suffix ("
                                                "fastq.gz or fastq). Please "
                                                "fix that issue." %
                                                input_path)

                # 3. Start fix_cutadapt
                fix_cutadapt = [self.get_tool('fix_cutadapt'),
                                temp_fifos["first_read_in"], 
                                temp_fifos["first_read_out"] ]
                if temp_fifos["second_read_in"] != None and \
                   temp_fifos["second_read_out"] != None:
                    fix_cutadapt.extend([
                        '--R2-in', temp_fifos["second_read_in"],
                        '--R2-out', temp_fifos["second_read_out"]
                    ])
                    
                exec_group.add_command(fix_cutadapt)
                
                # 4. Read data from first_read fifo
                with exec_group.add_pipeline() as fr_pigz_pipe:
                    # 4.1  command: Read from first_read fifos
                    cat = [self.get_tool('cat'),
                           temp_fifos["first_read_out"]]
                    # 4.2 Gzip output file
                    pigz = [self.get_tool('pigz'),
                            '--blocksize', '4096', 
                            '--processes', '2', 
                            '--stdout']
                    # 4.3 command: Write to output file in 4MB chunks
                    fr_stdout_path = run.add_output_file(
                        "%s" % read,
                        "%s%s.fastq.gz" %
                        (run_id, read_types["first_read"]),
                        input_paths)
                    dd = [self.get_tool('dd'),
                          'obs=4M',
                          'of=%s' % fr_stdout_path]
                            
                    fr_pigz_pipe.add_command(cat)
                    fr_pigz_pipe.add_command(pigz)
                    fr_pigz_pipe.add_command(dd)

                # 5. Read data from second_read fifo if there is one
                if temp_fifos["second_read_in"] != None and \
                   temp_fifos["second_read_out"] != None:
                    with exec_group.add_pipeline() as sr_pigz_pipe:
                        # 5.1  command: Read from first_read fifos
                        cat = [self.get_tool('cat'),
                               temp_fifos["second_read_out"]]
                        # 4.2 Gzip output file
                        pigz = [self.get_tool('pigz'),
                                '--blocksize', '4096', 
                                '--processes', '2', 
                                '--stdout']
                        # 4.3 command: Write to output file in 4MB chunks
                        sr_stdout_path = run.add_output_file(
                            "%s" % read,
                            "%s%s.fastq.gz" %
                            (run_id, read_types["second_read"]),
                            input_paths)
                        dd = [self.get_tool('dd'),
                              'obs=4M',
                              'of=%s' % sr_stdout_path]
                            
                        sr_pigz_pipe.add_command(cat)
                        sr_pigz_pipe.add_command(pigz)
                        sr_pigz_pipe.add_command(dd)
                    
#    def declare_runs(self):
#        # fetch all incoming run IDs which produce reads...
#        found_files = dict()
#        read_types = {'first_read': '-R1', 'second_read': '-R2'}
#        paired_end_info = dict()
#
#        for read in read_types.keys():
#            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/%s' % read):
#                if input_paths != [None]:
#                    paired_end_info[run_id] = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
#
#                # save information per file in found_files
#                if not run_id in found_files:
#                    found_files[run_id] = dict()
#
#                if not read in found_files[run_id]:
#                    found_files[run_id][read] = list()
#                # Check if we get exactly one input file
#                if len(input_paths) != 1:
#                    raise StandardError("Expected one input file.")
#                found_files[run_id][read].extend(input_paths)
#
#        for run_id in found_files.keys():
#            with self.declare_run(run_id) as run:
#                run.new_exec_group()
#                run.add_private_info('paired_end', paired_end_info[run_id])
#                for read in found_files[run_id].keys():
#                    run.add_output_file(read, 
#                        "%s-fixed%s.fastq.gz" % (run_id, read_types[read]), 
#                        found_files[run_id][read])
#            
#    def execute(self, run_id, run):
#        read_types = {'first_read': '-R1', 'second_read': '-R2'}
#        is_paired_end = run.get_private_info('paired_end')
#
#        with process_pool.ProcessPool(self) as pool:
#            first_read_path = run.get_single_output_file_for_annotation(
#                'first_read')
#            second_read_path = None
#            fifo_in_R1 = pool.get_temporary_fifo('fifo_in_R1', 'input')
#            fifo_in_R2 = None
#            fifo_out_R1 = pool.get_temporary_fifo('fifo_out_R1', 'output')
#            fifo_out_R2 = None
#            
#            with pool.Pipeline(pool) as pipeline:
#                cat = [self.get_tool('cat')]
#                cat.extend(
#                    run.get_input_files_for_output_file(first_read_path)
#                )
#                pigz = [self.get_tool('pigz'), '--decompress', '--processes', 
#                        '1', '--stdout']
#                pipeline.append(cat)
#                pipeline.append(pigz, stdout_path = fifo_in_R1)
#                
#            if is_paired_end:    
#                second_read_path = run.get_single_output_file_for_annotation(
#                    'second_read')
#                fifo_in_R2 = pool.get_temporary_fifo('fifo_in_R2', 'input')
#                fifo_out_R2 = pool.get_temporary_fifo('fifo_out_R2', 'output')
#                with pool.Pipeline(pool) as pipeline:
#                    cat = [self.get_tool('cat')]
#                    cat.extend(
#                        run.get_input_files_for_output_file(second_read_path)
#                        )
#                    pigz = [self.get_tool('pigz'), '--decompress', 
#                            '--processes', '1', '--stdout']
#                    pipeline.append(cat)
#                    pipeline.append(pigz, stdout_path = fifo_in_R2)
#        
#            fix_cutadapt = [self.get_tool('fix_cutadapt'), fifo_in_R1, 
#                            fifo_out_R1]
#            if is_paired_end:
#                fix_cutadapt.extend(['--R2-in', fifo_in_R2,
#                                     '--R2-out', fifo_out_R2])
#            pool.launch(fix_cutadapt)
#            
#            with pool.Pipeline(pool) as pipeline:
#                cat = [self.get_tool('cat'), fifo_out_R1]
#                pigz = [self.get_tool('pigz'), 
#                        '--blocksize', '4096', 
#                        '--processes', '2', 
#                        '--stdout']
#                pipeline.append(cat)
#                pipeline.append(pigz, stdout_path = first_read_path)
#                
#            if is_paired_end:
#                with pool.Pipeline(pool) as pipeline:
#                    cat = [self.get_tool('cat'), fifo_out_R2]
#                    pigz = [self.get_tool('pigz'), 
#                            '--blocksize', '4096', 
#                            '--processes', '2', 
#                            '--stdout']
#                    pipeline.append(cat)
#                    pipeline.append(pigz, stdout_path = second_read_path)
