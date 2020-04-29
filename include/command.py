from uaperrors import UAPError
import os
from logging import getLogger
import pipeline_info
import exec_group

logger = getLogger('uap_logger')


def abs2rel_path(func):
    '''
    A decoraror function to replace absolute paths with relative paths.
    It also removes the deprectaed output path placeholders for
    backwards compatibility with old step implementation.
    '''

    def inner(self, *args):
        run = self.get_run()
        working_dir = run.get_temp_output_directory()
        abs_dest = run.get_step().get_pipeline().config['destination_path']
        rel_path = os.path.relpath(abs_dest, working_dir)

        def repl(text):
            if isinstance(text, str):
                if text.endswith(abs_dest):
                    return text.replace(abs_dest, rel_path)
                return text.replace(abs_dest+os.sep, rel_path+os.sep)
            elif isinstance(text, list) or isinstance(text, set):
                return [repl(element) for element in text]
            elif text is None:
                return None
            else:
                raise UAPError("Function %s does not return string or "
                               "list of strings." % func.__name__)
        return repl(func(self, *args))
    return inner


def quote(cmd_args):
    """
    An arument list is combined into a string and arguments
    containing bash special characters are single quoted.
    """
    if not isinstance(cmd_args, str):
        return ' '.join(quote(c) for c in cmd_args)
    if "'" in cmd_args:
        cmd_args = "'%s'" % cmd_args.replace("'", "\\'")
    elif any(s in cmd_args for s in ' |*+";?&()[]<>$#`\t\n'):
        cmd_args = "'%s'" % cmd_args
    return cmd_args


class CommandInfo(object):
    def __init__(self, eop, command, stdout_path=None, stderr_path=None):
        # eop = exec_group or pipeline
        self._eop = eop
        self._command = list()
        self._stdout_path = stdout_path
        self._stderr_path = stderr_path
        self._tool = str
        self._output_files_per_connection = dict()

        for _ in command:
            if _ == command[0]:
                self._tool = _
                if isinstance(command[0], list):
                    self._tool = _[-1]
                pass
            elif not isinstance(_, str):
                raise TypeError(
                    "Non-string element %s in command %s" %
                    (_, command))
            self._command.append(_)

    def get_run(self):
        if isinstance(self._eop, pipeline_info.PipelineInfo):
            run_info = self._eop.get_exec_group().get_run()
        elif isinstance(self._eop, exec_group.ExecGroup):
            run_info = self._eop.get_run()
        else:
            run_info = None
        return run_info

    def set_command(self, command):
        if not isinstance(command, list):
            raise UAPError("Given non-list command: %s" % command)
        self._command = command

    def set_stdout_path(self, stdout_path):
        if stdout_path is not None:
            self._stdout_path = stdout_path

    @abs2rel_path
    def get_stdout_path(self):
        return self._stdout_path

    def set_stderr_path(self, stderr_path):
        if stderr_path is not None:
            self._stderr_path = stderr_path

    @abs2rel_path
    def get_stderr_path(self):
        return self._stderr_path

    @abs2rel_path
    def get_command(self):
        return self._command

    def get_command_string(self, replace_path=False):
        '''
        Returns a string representation of the command that is runable in bash.
        '''
        # We cannot use subprocess.list2cmdline here since it does not
        # quote for bash and leaves special characters such as | as is.
        out = ''
        if self._stdout_path:
            out += ' > %s' % self._stdout_path
        if self._stderr_path:
            out += ' 2> %s' % self._stderr_path
        if replace_path is not True:
            return quote(self.get_command()) + out
        # replace tool call with its name
        cmd = self.get_command()
        map = self.get_run().get_step().get_path_tool()
        tool = map.get(' '.join(cmd[0]))
        if tool is None:
            tool = cmd[0]
            for path, tool in map.items():
                tool = tool.replace(path, tool)
        return quote([tool] + cmd[1:]) + out
