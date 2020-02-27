from uaperrors import UAPError
import sys
import os
from logging import getLogger
import pipeline_info
import exec_group

logger=getLogger('uap_logger')

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
                raise TypeError("Non-string element %s in command %s" % (_, command))
            self._command.append(_)

    def replace_output_dir_du_jour(func):
        def inner(self, *args):
            run = self.get_run()
            # Collect info to remove derecated placeholders
            # and replace absolute paths with relative paths
            placeholder = run.get_output_directory_du_jour_placeholder()
            placeholder += '' if placeholder.endswith('/') else '/'
            pl_slashless = placeholder[:-1]

            working_dir = run.get_temp_output_directory()
            abs_dest = run.get_step().get_pipeline().config['destination_path']
            rel_path = os.path.relpath(abs_dest, working_dir)
            def repl(text):
                if isinstance(text, str):
                    text = text.replace(abs_dest, rel_path)
                    text = text.replace(placeholder, '')
                    text = text.replace(pl_slashless, '.')
                    return text
                elif isinstance(text, list) or isinstance(text, set):
                    return [repl(element) for element in text]
                elif text is None:
                    return None
                else:
                    raise UAPError("Function %s does not return string or "
                                   "list of strings." % func.__name__)
            return repl(func(self, *args))
        return inner

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
        if stdout_path != None:
            self._stdout_path = stdout_path

    @replace_output_dir_du_jour
    def get_stdout_path(self):
        return self._stdout_path

    def set_stderr_path(self, stderr_path):
        if stderr_path != None:
            self._stderr_path = stderr_path

    @replace_output_dir_du_jour
    def get_stderr_path(self):
        return self._stderr_path

    @replace_output_dir_du_jour
    def get_command(self):
        return self._command
