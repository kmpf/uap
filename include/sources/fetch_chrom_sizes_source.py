from uaperrors import StepError
import sys
import logging
import traceback
import os

from abstract_step import *

logger = logging.getLogger("uap_logger")


class FetchChromSizes(AbstractStep):

    def __init__(self, pipeline):
        super(FetchChromSizesSource, self).__init__(pipeline)

        self.add_connection('out/chromosome_sizes')

        self.require_tool('fetchChromSizes')
        self.require_tool('cp')

        self.add_option('ucsc-database', str, optional=False,
                        description="Name of UCSC database e.g. hg38, mm9")
#         self.add_option('path', str, optional=False,
#                         description="directory to move file to")

    def runs(self, run_ids_connections_files):
        '''

        '''
        output_dir = self.get_option('path')
        if not os.path.exists(output_dir):
            exc_type, exc_value, exc_traceback = sys.exc_info()

            tb = "Stack trace:"
            for stack_entries in traceback.extract_tb(exc_traceback, 4):
                tb += "\n"
                tb += " %s, line %s, %s %s" % (stack_entries)

            logger.debug(tb)
            raise StepError(self,
                'Output directory (%s) does not exist. Please create it.' %
                output_dir)

        ucsc_database = self.get_option('ucsc-database')

        output_filename = "%s.chrom.sizes" % ucsc_database

#         out_file = os.path.abspath(os.path.join(output_dir, output_filename))

        # Declare a new run
        with self.declare_run(ucsc_database) as run:
            temp_filename = run.add_temporary_file(suffix=output_filename)
            with run.new_exec_group() as exec_group:
                fetch_chrom_sizes = [
                    self.get_tool('fetchChromSizes'),
                    ucsc_database
                ]
                exec_group.add_command(
                    fetch_chrom_sizes,
                    stdout_path=temp_filename
                )

                cp = [self.get_tool('cp'), '--update', temp_filename, output_filename]
                exec_group.add_command(cp)

                run.add_output_file('chromosome_sizes', output_filename, [])
