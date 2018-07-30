import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SamtoolsMerge(AbstractStep):

    '''The step samtools_merge builds on 'samtools merge' to merge sorted SAM/BAM files and output
    SAM/BAM/CRAM files.

    Merge adds readgroup tag to preserve the information on the orginial sample.

    This step wraps samtools merge from release 1.7.

    Documentation:

        http://www.htslib.org/doc/samtools.html

    '''

    def __init__(self, pipeline):
        super(SamtoolsMerge, self).__init__(pipeline)

        self.set_cores(6)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log')
        self.add_connection('out/err')

        self.require_tool('dd')
        self.require_tool('samtools')
        self.require_tool('pigz')

        # [Options for 'samtools':]
        self.add_option('n', bool, optional=True,
                        description='sort by read names')
        self.add_option('t', str, optional=True,
                        description='Input files are sorted by TAG value')
        self.add_option('r', bool, optional=True,
                        description='attach RG tag (inferred from file names)')
        self.add_option('u', bool, optional=True,
                        description='uncompressed BAM output')
        self.add_option('f', bool, optional=True,
                        description='overwrite the output BAM if exist')
        self.add_option('1', bool, optional=True,
                        description='compress level 1')
        self.add_option('l', int, optional=True,
                        description='compression level, from 0 to 9 [-1]')
        self.add_option('threads', int, optional=True,
                        description='Number of additional threads to use [0]')
        self.add_option('R', str, optional=True,
                        description='merge file in the specified region STR [all]')
        self.add_option('c', bool, optional=True,
                        description='Combine @RG headers with colliding IDs [alter IDs to be distinct]')
        self.add_option('p', bool, optional=True,
                        description='Combine @PG headers with colliding IDs [alter IDs to be distinct]')
        self.add_option('s', str, optional=True,
                        description='override random seed')
        # options -h (which header is copied to output) and -b a list with input bam filenames are
        # not implemented, because it is taken care of inside uap
        # in/output format options are not implemented. We stick to BAM for now.

        self.add_option('run_id', str, optional=True, default="mergeSBam",
                        description="A name for the run. Since this step merges multiple samples "
                        "into a single one, the run_id cannot be the sample name anymore.")

        # [Options for 'dd' and 'pigz':]
        self.add_option('dd-blocksize', str, optional = True, default = "4M")
        self.add_option('pigz-blocksize', str, optional = True, default = "4096")

    def runs(self, run_ids_connections_files):

        options = ['n','r','u','f','1','l','R','c','p','s']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()

        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))

        if self.is_option_set_in_config('threads'):
            option_list.append('-@ %d' % self.get_option('threads'))

        run_id = self.get_option('run_id')

        # Get all input .bam files that should be merged
        input_list = list()

        for sample_id in run_ids_connections_files.keys():

            input_paths = run_ids_connections_files[sample_id]["in/alignments"]
            input_list.append(input_paths[0])


        with self.declare_run(run_id) as run:

            out_file = run.add_output_file('alignments',
                                           '%s-samtools-merged.bam' % run_id,
                                           input_paths)
            # log files
            log_outfile = run.add_output_file('log',
                                              '%s-samtools-merged-stdout.txt' % run_id,
                                              input_paths)
            ## let's see what's provided here
            err_outfile = run.add_output_file('err',
                                              '%s-samtools-merged-stderr.txt' % run_id,
                                              input_paths)

            with run.new_exec_group() as exec_group:

                samtools_merge = [self.get_tool('samtools'), 'merge']
                samtools_merge.extend(option_list)
                samtools_merge.extend([out_file])
                samtools_merge.extend(input_list)

                exec_group.add_command(samtools_merge,
                                       stdout_path = log_outfile,
                                       stderr_path = err_outfile)

