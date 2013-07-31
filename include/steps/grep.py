import sys
from abstract_step import *
import pipeline
import io_step

class Grep(io_step.IOStep):
    
    def __init__(self, pipeline):
        super(Grep, self).__init__(pipeline, 'grep')
        self.add_option('pattern', str, optional=False)

    def tool_command_line(self):
        return [self.option('pattern')]