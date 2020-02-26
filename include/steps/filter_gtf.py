import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class filterGTF(AbstractStep):

    '''
    custom script to filter merged genocde gtf from cufflinks or stringtie by classcode
    '''


    def __init__(self, pipeline):
        super(filterGTF, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/assembling')
        self.add_connection('out/assembling')
        self.add_connection('out/log_stderr')
        self.require_tool('filter_gtf')

        self.add_option('remove-unstranded', bool, optional=False,
                        description='Removes transcripts without strand specifity')

        self.add_option('class-code-only-in-transcript-feature', bool, optional=False,
                        description='Removes transcripts without strand specifity')

        self.add_option('remove-unwanted-chr', bool, optional=False,
                        description='keeps chr1 ..2 and chrX, chrY, chrMT')

        self.add_option('string', str, optional=True,
                        description='string to match in gtf field gene_name for discarding')

        self.add_option('remove-by-field-match', str, optional=True,
                        description='select gft field like gene_id, gene_name which will match against --string')

        self.add_option('keep-by-class', bool, optional=False,
                        description='"keep  gtf if any class is found in class_code field, requieres class-list-keep')
        self.add_option('class-list-keep', str, optional=True,
                        description="class codes to be kept possible '=,c,j,e,i,o,p,r,u,x,s,.'")



    def runs(self, run_ids_connections_files):

        options=['remove-unstranded', 'class-code-only-in-transcript-feature',
                 'remove-unwanted-chr', 'string', 'remove-by-field-match',
                 'keep-by-class', 'class-list-keep']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append( '--%s' % option )
                option_list.append( str(self.get_option(option)) )


        for run_id in run_ids_connections_files.keys():

             with self.declare_run(run_id) as run:
                input_path  = run_ids_connections_files[run_id]['in/assembling'][0]


                assembling = run.add_output_file(
                    'assembling',
                    '%s.annotated.gtf' % run_id,
                    [input_path])

                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-log_stderr.txt' % run_id,
                    [input_path])


                # check, if only a single input file is provided
                len_input = run_ids_connections_files[run_id]['in/assembling']
                if len(len_input) != 1:
                    raise StandardError("Expected exactly one assembling file, but got this: %s" % len_input)

                with run.new_exec_group() as exec_group:

                   # exec_group.add_command(loci, stdout_path = loci_file)

                    with exec_group.add_pipeline() as pipe:
                        filter_gtf = [self.get_tool('filter_gtf')]
                        filter_gtf.extend(option_list)
                        filter_gtf.append(input_path)

                        pipe.add_command(filter_gtf, stdout_path=assembling, stderr_path=log_stderr)

