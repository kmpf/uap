"""

"""
import command as command_info


class PipelineInfo(object):
    def __init__(self, exec_group):
        self._exec_group = exec_group
        self._commands = list()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def add_command(self, command, stdout_path=None, stderr_path=None):
        command = command_info.CommandInfo(
            self, command, stdout_path=stdout_path,
            stderr_path=stderr_path
        )
        self._commands.append(command)
        return command

    def get_commands(self):
        return self._commands

    def get_command_string(self, replace_path=False):
        return ' | '.join(c.get_command_string(replace_path=replace_path)
                          for c in self.get_commands())

    def get_exec_group(self):
        return self._exec_group
