"""

"""
import os

import abstract_step
import command as command_info
import misc
import pipeline_info
import run

class ExecGroup(object):
    def __init__(self, run):
        self._run = run
        self._pipes_and_commands = list()
        self._fifos = dict()

    def __enter__(self):        
        return self

    def __exit__(self, type, value, traceback):
        pass
        
    def new_pipeline(self):
        pipeline = pipeline_info.PipelineInfo(self)
        self._pipes_and_commands.append(pipeline)
        return pipeline

    def new_command(self, command, stdout_path=None, stderr_path=None):
        command = command_info.CommandInfo(self, command, 
                                           stdout_path = stdout_path, 
                                           stderr_path = stderr_path
                                       )
        self._pipes_and_commands.append(command)
        return command

    def get_pipes_and_commands(self):
        return self._pipes_and_commands

    def get_run(self):
        return self._run

    def new_temporary_fifo(self, prefix, designation = None):
        placeholder = "<temp-fifo-%s-%s>" % (prefix, designation)
        self._fifos[placeholder] = {'prefix' : prefix,
                                    'designation': designation}
        return placeholder

    def get_temporary_fifo(self, placeholder):
        return self._fifos[placeholder]
