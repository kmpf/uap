from uaperrors import UAPError
import logging
from abstract_step import AbstractStep
import os
import sys

logger = logging.getLogger('uap_logger')


class FastxReverseComplement(AbstractStep):
    '''
    wrapper class for fastx_reverse_complement from fastx toolkit
    creates reverse complement of fasta and fastq files.
    http://hannonlab.cshl.edu/fastx_toolkit/
    '''

    def __init__(self, pipeline):
        super(FastxReverseComplement, self).__init__(pipeline)

        self.set_cores(1)  # muss auch in den Decorator

        self.add_connection('in/fastx')
        self.add_connection('out/fastx')

        self.require_tool('fastx_reverse_complement')
        self.require_tool('pigz')
        self.require_tool('cat')
    
        self.add_option('prefix', str, default=None, optional=True,
                        description="Add Prefix to sample name")

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            new_run_id = run_id
            if self.is_option_set_in_config('prefix'):
               new_run_id = self.get_option('prefix') + '_' + run_id

            with self.declare_run(new_run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/fastx']
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise UAPError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                out = run.add_output_file(
                    "fastx",
                    "%s_%s-fastq.gz" %  (new_run_id, 'R1'),
                    input_paths) 
       
                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
          
                    # 1.1 command: Uncompress file
                        if is_gzipped:
                            pigz = [self.get_tool('pigz'),
                                    '--decompress',
                                    '--processes', '1',
                                    '--stdout']
                            pigz.extend(input_paths)
                            pipe.add_command(pigz)
                           
                        else:
                            cat = [self.get_tool('cat')]
                            cat.extend(input_paths)
                            pipe.add_command(cat)
                                   

                        # 1. Run  fastx  for input file
                        fastx_revcom = [self.get_tool('fastx_reverse_complement')]
                        # gzip 
                        fastx_revcom.extend(['-z'])
                        pipe.add_command(fastx_revcom,  stdout_path=out )
                        

