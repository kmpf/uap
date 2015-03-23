"""

"""


class PipelineInfo(object):
    def __init__(self, exec_group):
        self._exec_group = exec_group
        self._commands = list()

    def new_command(self):
        command = Command_info(self)
        self._commands.append(command)
        return command
        
    def get_commands(self):
        return self._commands

    def get_exec_group(self):
        return self._exec_group
