from uaperrors import UAPError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Samtools(AbstractStep):
    '''
    The step samtools wraps parts of the 'samtools' packages. It is intended for 
    reformatting SAM/BAM files and not completely implemented.

    Feel free to add the samtools options you need!

    The options listed below are implemented.
    
    For a description/explanation about SAM flags, we refer to 
    - the samtools manual page http://www.htslib.org/doc/samtools.html
    - the Picard page to explain SAM flags https://broadinstitute.github.io/picard/explain-flags.html

    '''

    def __init__(self, pipeline):
        super(Samtools, self).__init__(pipeline)
        
        self.set_cores(8)
        
        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        
        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('pigz')


        
        # general samtools options
        # -h 
        self.add_option('keep_header', bool, optional = True, default = True,
                        description = 'Include the header in the output.') 
        # -b
        self.add_option('output_bam', bool, optional = True, default = False,
                        description = 'Output in the BAM format.')
        # -t <FILE>
        self.add_option('genome-faidx', str, optional = False)

        # Quality 
        # -q <INT>
        self.add_option('q_mapq', int, optional = True, default = 0,
                        description = 'Skip alignments with MAPQ smaller than '
                        'this value.')

        ### samtools view
        # Filter 
        # -f <INT>
        self.add_option('f_keep', int, optional = True, default = 0,
                        description = 'Only output alignments with all bits set in '
                        'INT present in the FLAG field.')
        # -F <INT>
        self.add_option('F_skip', int, optional = True, default = 0,
                        description = 'Do not output alignments with any bits set '
                        'in INT present in the FLAG field.')

        ### samtools sort
        # Sorting        
        # -n 
        self.add_option('sort-by-name', bool, default = False,
                        description = 'Sort by read names (i.e., the QNAME field) '
                        'rather than by chromosomal coordinates.')
        self.add_option('temp-sort-dir', str, optional = False,
                        description = 'Intermediate sort files are stored into'
                        'this directory.')

        # [Options for 'dd':]
        self.add_option('dd-blocksize', str, optional = True, default = "2048",
                        description = 'Blocksize for dd tool.')

    def runs(self, run_ids_connections_files):
        
        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise UAPError("Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                                 in ['.gz', '.gzip'] else False

                if self.is_option_set_in_config('temp-sort-dir'):
                    if not os.path.isdir(self.get_option('temp-sort-dir')):
                        #dir not present
                        raise UAPError("Directory %s not found" % self.get_option('temp-sort-dir'))
                    if not os.access(self.get_option('temp-sort-dir'), os.W_OK):
                        #not accessible
                        raise UAPError("Directory %s not accessible." % self.get_option('temp-sort-dir'))
                

                    with run.new_exec_group() as exec_group:

                        with exec_group.add_pipeline() as pipe:
                            # 1. command: Read file in chunks
                            dd_in = [
                                self.get_tool('dd'),
                                'ibs=%s' % self.get_option('dd-blocksize'),
                                'if=%s' % input_paths[0]
                            ]
                            pipe.add_command(dd_in)

                            # 1.1 command: Uncompress file to fifo
                            if is_gzipped:
                                pigz = [self.get_tool('pigz'),
                                        '--decompress',
                                        '--processes', str(self.get_cores()),
                                        '--blocksize', self.get_option('dd-blocksize'),
                                        '--stdout']
                                pipe.add_command(pigz)

                            # 2. command: Convert sam to bam
                            samtools_view = [
                                self.get_tool('samtools'), 'view']
                            
                            # 3. command: Sort BAM input
                            samtools_sort = [
                                self.get_tool('samtools'), 'sort']

                            if self.get_option('sort-by-name'):
                                samtools_sort.append('-n')
                            
                            # add the general options
                            # header
                            if self.get_option('keep_header'):
                                samtools_view.append('-h')

                            # output
                            if self.get_option('output_bam'):
                                samtools_view.append('-b')
                                samtools_sort.extend(['-O', 'bam'])
                            # quality
                            if self.get_option('q_mapq') > 0:
                                samtools_view.extend(['-q', 
                                                     str(self.get_option('q_mapq'))])
                                # no need to do this for samtools_sort

                            # filter options
                            # keep
                            if self.get_option('f_keep') > 0:
                                samtools_view.extend(['-f',
                                                      str(self.get_option('f_keep'))])
                            # skip
                            if self.get_option('F_skip') > 0:
                                samtools_view.extend(['-F',
                                                      str(self.get_option('F_skip'))])

                            samtools_view.extend(['-t', self.get_option('genome-faidx'), 
                                                  '-', '-@', str(self.get_cores())])
                            samtools_sort.extend(['-T', 
                                                  os.path.join(self.get_option('temp-sort-dir'), run_id), 
                                                  '-', '-@', str(self.get_cores())]
                            )
                            
                            # add to pipe
                            pipe.add_command(samtools_view)

                            pipe.add_command(samtools_sort)

                            # 4. command: dd
                            dd_out = [
                                self.get_tool('dd'),
                                'obs=%s' % self.get_option('dd-blocksize')
                            ]
                            pipe.add_command(
                                dd_out,
                                stdout_path = run.add_output_file(
                                    'alignments',
                                    '%s.samtools.bam' % run_id,
                                    input_paths)
                            )
