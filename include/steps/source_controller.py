from uaperrors import UAPError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class SourceController(AbstractStep):
    '''
    This step combines all inputs, produces a symlink to each file
    and hashes them. It may be use to inherit from source steps
    so changes in the source files can be detected later on.
    '''

    def __init__(self, pipeline):
        super(SourceController, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/raw')
        self.add_connection('out/raw')

        self.require_tool('ln')

        self.add_option('cores', int, optional=True, default=4,
                description='Number of threads used to calculate the hash sums.')

    def runs(self, cc):

        group = self.declare_run('links')
        execg = group.new_exec_group()
        for run_id, files in cc.connection_items('in/raw'):
            for file in files:
                link_name = run_id + '-' + os.path.basename(file)
                link = group.add_output_file('raw', link_name, [file])
                ln = [self.get_tool('ln'), '-s', file, link_name]
                execg.add_command(ln)

