"""

"""
import os

import abstract_step
from command import CommandInfo
import misc
from pipeline_info import PipelineInfo
import run

class ExecGroup(object):
    def __init__(self, run):
        self._run = run
        self._pipes_and_commands = list()
        
    def new_pipeline(self):
        pipeline = pipeline.PipelineInfo(self)
        self._pipes_and_commands.append(pipeline)
        return pipeline

    def new_command(self, command, stderr_path=None, stdout_path=None):
        command = CommandInfo(self, command, stderr_path,
                                            stdout_path)
        self._pipes_and_commands.append(command)
        return command

    def get_pipes_and_commands(self):
        return self._pipes_and_commands

    def get_run(self):
        return self._run
