import sys
from abstract_step import *
import pipeline

class Break(AbstractStep):
    '''
    This steps breaks the processing flow and does nothing. Because it
    returns no runs, its children never get run. Use it to conveniently
    cut off branches of the step tree (it's like an uncommenting feature).
    '''

    def __init__(self, pipeline):
        super(Break, self).__init__(pipeline)

    def setup_runs(self, input_run_info):
        return {}

    def execute(self, run_id, run_info):
        pass
