import sys
from abstract_step import *
import pipeline

'''
This steps breaks the processing flow and does nothing. Because it
returns no runs, its children never get run.
'''

class Break(AbstractStep):
    def __init__(self, pipeline):
        super(Break, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        return {}

    def execute(self, run_id, run_info):
        pass
