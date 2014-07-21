import sys
from abstract_step import *
import pipeline
import io_step

class Head(io_step.IOStep):
    
    def __init__(self, pipeline):
        super(Head, self).__init__(pipeline, 'head')
        self.add_option('lines', int, default = 1000)

    def tool_command_line(self):
        return ['-n', str(self.get_option('lines'))]
