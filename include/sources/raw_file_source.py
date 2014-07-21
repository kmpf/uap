import sys
from abstract_step import *
import os

class RawFileSource(AbstractSourceStep):

    def __init__(self, pipeline):
        super(RawFileSource, self).__init__(pipeline)
        self.add_option('path', str, optional=False)
        self.add_connection('out/raw')

    def declare_runs(self):
        with self.declare_run('nop') as run:
            run.add_output_file('raw', os.path.abspath(self.get_option('path')), [])
