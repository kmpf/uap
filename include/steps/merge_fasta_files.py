import os
import sys
from logging import getLogger
import os
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class MergeFastaFiles(AbstractStep):
    '''
    This step merges all .fasta(.gz) files belonging to a certain sample.
    The output files are gzipped.
    '''
    
    def __init__(self, pipeline):
        super(MergeFastaFiles, self).__init__(pipeline)
        
        self.set_cores(12)
        
        self.add_connection('in/sequence')
        self.add_connection('out/sequence')

        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('mkfifo')
        self.require_tool('pigz')

        self.add_option('compress-output', bool, optional=True, default=True,
                        description="If set to true output is gzipped.")
        self.add_option('merge-all-runs', bool, optional=True,
                        default=False, description="If set to true sequences "
                        "from all runs are merged")
        self.add_option('output-fasta-basename', str, optional=True, default="",
                        description="Name used as prefix for FASTA output.")

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "256k")

    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        run_ids = set(run_ids_connections_files.keys())
        for run_id in run_ids_connections_files.keys():
            input_paths = list()

            # Do everything that's necessary for 'merge-all-runs'
            if self.get_option('merge-all-runs'):
                run_ids.remove(run_id)
                input_paths.extend(
                    run_ids_connections_files[run_id]['in/sequence'])
                if len(run_ids) > 0: continue
                run_id = 'all_sequences'
            else:
                input_paths = run_ids_connections_files[run_id]['in/sequence']

            fasta_basename = run_id
            if self.get_option('output-fasta-basename'):
                fasta_basename = "%s-%s" % (
                    self.get_option('output-fasta-basename'), run_id)

            with self.declare_run(fasta_basename) as run:

                if input_paths == [None]:
                    run.add_empty_output_connection("sequence")
                else:
                    temp_fifos = list()
                    exec_group = run.new_exec_group()
                    for input_path in input_paths:
                        # Gzipped files are unpacked first
                        # !!! Might be worth a try to use fifos instead of
                        #     temp files!!!
                        # 1. Create temporary fifo
                        temp_fifo = run.add_temporary_file(
                            "fifo-%s" %
                            os.path.basename(input_path) )
                        temp_fifos.append(temp_fifo)
                        mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                        exec_group.add_command(mkfifo)

                        is_gzipped = True if os.path.splitext(input_path)[1]\
                                     in ['.gz', '.gzip'] else False

                        # 2. Output files to fifo
                        if is_gzipped:
                            with exec_group.add_pipeline() as unzip_pipe:
                                # 2.1 command: Read file in 4MB chunks
                                dd_in = [
                                    self.get_tool('dd'),
                                    'ibs=%s' % self.get_option('dd-blocksize'),
                                    'if=%s' % input_path
                                ]
                                unzip_pipe.add_command(dd_in)

                                # 2.2 command: Uncompress file to fifo
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--stdout']
                                unzip_pipe.add_command(pigz)

                                # 2.3 Write file in 4MB chunks to fifo
                                dd_out = [
                                    self.get_tool('dd'),
                                    'obs=%s' % self.get_option('dd-blocksize'),
                                    'of=%s' % temp_fifo
                                ]
                                unzip_pipe.add_command(dd_out)
                        
                        elif os.path.splitext(input_path)[1] in\
                             ['.fastq', '.fq', '.fasta', '.fa', '.fna']:
                            # 2.1 command: Read file in 4MB chunks and
                            #              write to fifo in 4MB chunks
                            dd_in = [
                                self.get_tool('dd'),
                                'bs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_path,
                                'of=%s' % temp_fifo
                            ]
                            exec_group.add_command(dd_in)
                        else:
                            logger.error("File %s does not end with any "
                                         "expected suffix (fastq.gz or fastq)."
                                         " Please fix that issue." %
                                         input_path)
                            sys.exit(1)
                    # 3. Read data from fifos
                    with exec_group.add_pipeline() as pigz_pipe:
                        # 3.1 command: Read from ALL fifos
                        cat = [self.get_tool('cat')]
                        cat.extend(temp_fifos)
                        pigz_pipe.add_command(cat)

                        # 3.2 Gzip output file
                        out_file = "%s.fasta" % fasta_basename
                        if self.get_option('compress-output'):
                            out_file = "%s.fasta.gz" % fasta_basename
                            pigz = [self.get_tool('pigz'),
                                    '--stdout']
                            pigz_pipe.add_command(pigz)

                        # 3.3 command: Write to output file in 4MB chunks
                        stdout_path = run.add_output_file(
                            "sequence",
                            out_file,
                            input_paths)
                        dd = [
                            self.get_tool('dd'),
                            'obs=%s' % self.get_option('dd-blocksize'),
                            'of=%s' % stdout_path
                        ]
                        pigz_pipe.add_command(dd)
