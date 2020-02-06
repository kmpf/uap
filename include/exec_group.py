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

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def add_pipeline(self):
        pipeline = pipeline_info.PipelineInfo(self)
        self._pipes_and_commands.append(pipeline)
        return pipeline

    def add_command(self, command, stdout_path=None, stderr_path=None):
        try:
            command = command_info.CommandInfo(self, command,
                                               stdout_path = stdout_path,
                                               stderr_path = stderr_path
                                           )
        except TypeError, e:
            raise UAPError('During declaration of step "%s": %s' %
                    (str(self._run.get_step()), e.message))
        self._pipes_and_commands.append(command)
        return command

    def get_pipes_and_commands(self):
        return self._pipes_and_commands

    def get_run(self):
        return self._run
