import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class Stringtie(AbstractStep):

    def __init__(self, pipeline):
        super(Stringtie, self).__init__(pipeline)

        self.set_cores(12)

        self.add_connection('in/alignments')
        self.add_connection('out/assembling')
        self.add_connection('out/gene_abund')
        self.add_connection('out/cov_refs')
        self.add_connection('out/log_stderr')

        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('stringtie')

        self.add_option('G', str, optional=False,
                        description="use reference transcript annotation to guide assembly")
        self.add_option('v', bool, optional=True,
                        description='Turns on verbose mode, printing bundle processing details')
        self.add_option('p', int, default=1, optional=True,
                        description='Specify the number of processing threads (CPUs) to use for transcript assembly. The default is 1')
        self.add_option('m', int, optional=True,
                        description='Sets the minimum length allowed for the predicted transcripts. Default: 200')
        self.add_option('l', str, optional=True,
                        description='Sets <label> as the prefix for the name of the output transcripts. Default: STRG')
        self.add_option('f', float, optional=True,
                        description='Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts. Default: 0.1')

        self.add_option('fr', bool, optional=True,
                    description='assume stranded library fr-secondstrand')

        self.add_option('rf', bool, optional=True,
                    description='assume stranded library fr-firststrand')


    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('p'))

        options=['G', 'v','p', 'm', 'l', 'f']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append( '-%s' % option )
                option_list.append( str(self.get_option(option)) )
        
        if (self.is_option_set_in_config('fr') and self.get_option('fr')):
            option_list.append('--fr')

        if (self.is_option_set_in_config('rf') and self.get_option('rf')):
            option_list.append('--rf')
            

        for run_id in run_ids_connections_files.keys():

             with self.declare_run(run_id) as run:
                alignments  = run_ids_connections_files[run_id]['in/alignments'][0]
                input_paths = alignments


                assembling = run.add_output_file(
                    'assembling',
                    '%s-assembling.gtf' % run_id,
                    [input_paths])

                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-log_stderr.txt' % run_id,
                    [input_paths])

                gene_abund = run.add_output_file(
                    'gene_abund',
                    '%s-gene_abund.tab' % run_id,
                    [input_paths])

                cov_refs = run.add_output_file(
                    'cov_refs',
                    '%s-cov_refs.gtf' % run_id,
                    [input_paths])

                # check reference annotation
                if not os.path.isfile(self.get_option('G')):
                    logger.error(
                        "The path %s provided to option 'G' is not a file."
                        % self.get_option('G') )
                    sys.exit(1)

                # check, if only a single input file is provided
                len_input = run_ids_connections_files[run_id]['in/alignments']
                if len(len_input) != 1:
                    raise StandardError("Expected exactly one alignments file., but got this %s" % input_paths)

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
                        stringtie = [self.get_tool('stringtie'), alignments, '-o', assembling,
                                     '-A', gene_abund, '-C', cov_refs
                        ]
                        stringtie.extend(option_list)
                        pipe.add_command(stringtie, stdout_path=assembling, stderr_path=log_stderr)

