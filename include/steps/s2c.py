from uaperrors import StepError
import sys
from abstract_step import *
import process_pool
import yaml
import os
from logging import getLogger

logger=getLogger('uap_logger')

class S2C(AbstractStep):
    '''
    s2c formats the output of segemehl mapping to be compatible with
    the cufflinks suite of tools for differential expr. analysis of
    RNA-Seq data and their visualisation.
    For details on cufflinks we refer to the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/

    '''
    def __init__(self, pipeline):
        super(S2C, self).__init__(pipeline)

        self.set_cores(6)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/log')

        self.require_tool('s2c')
        self.require_tool('fix_s2c')
        self.require_tool('samtools')
        self.require_tool('pigz')
        self.require_tool('cat')
        self.require_tool('dd')

        self.add_option('tmp_dir', str, optional=False,
                        description="Temp directory for 's2c.py'. This can be "
                        "in the /work/username/ path, since it is only "
                        "temporary.")
        self.add_option('maxDist', int, optional=True,
                        description="specifies the maximal distance of a splice junction. "
                        "junctions with disctance higher than this value are classified as "
                        "fusions (default is 200.000nt)")

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/alignments']
                # check, if only a single input file is provided
                if len(input_paths) != 1:
                    raise Exception("Expected exactly one alignments file., but got this %s" % input_paths)

                if self.is_option_set_in_config('tmp_dir'):
                    if not os.path.isdir(self.get_option('tmp_dir')):
                        #dir not present
                        raise StepError(self, "Directory %s not found" % self.get_option('tmp_dir'))
                    if not os.access(self.get_option('tmp_dir'), os.W_OK):
                        #not accessible
                        raise StepError(self, "Directory %s not accessible." % self.get_option('tmp_dir'))

                alignments_path = input_paths[0]
                cat = [self.get_tool('cat'), alignments_path]
#                pigz = [self.get_tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                pigz = [self.get_tool('pigz'), '--decompress', '--processes', str(self.get_cores()), '--stdout']
                s2c = [self.get_tool('s2c'), '-s', '/dev/stdin', '-o', self.get_option('tmp_dir')]
                if self.is_option_set_in_config('maxDist'):
                    s2c.extend(['-d', str(self.get_option('maxDist'))])

                fix_s2c = [self.get_tool('fix_s2c')] # schreibt .sam nach stdout
#                pigz2 = [self.get_tool('pigz'), '--processes', '2', '--stdout']
                pigz2 = [self.get_tool('pigz'), '--processes', str(self.get_cores()), '--stdout']

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as s2c_pipe:
                        s2c_pipe.add_command(cat)
                        s2c_pipe.add_command(pigz)
                        s2c_pipe.add_command(s2c)
                        s2c_pipe.add_command(fix_s2c)
                        s2c_pipe.add_command(pigz2,
                                             stdout_path= run.add_output_file(
                                                 'alignments',
                                                 '%s-cufflinks-compatible.sam.gz' % run_id,
                                                 input_paths))
