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
                logger.error("Non-string element %s in command %s" % (_, command))
                sys.exit(1)
            self._command.append(_)

    def replace_output_dir_du_jour(func):
        def inner(self, *args):
            run_info = self.get_run()
            # Collect info to replace du_jour placeholder with temp_out_dir
            placeholder = run_info.get_output_directory_du_jour_placeholder()
            temp_out_dir = run_info.get_output_directory_du_jour()

            command = None
            ret_value = func(self, *args)
            if ret_value is None:
                return(None)
            if isinstance(ret_value, list) or isinstance(ret_value, set):
                command = list()
                for string in ret_value:
                    if string != None and placeholder in string and\
                       isinstance(temp_out_dir, str):
                        command.append(
                            string.replace(placeholder, temp_out_dir))
                    else:
                        command.append(string)
            elif isinstance(ret_value, str):
                if ret_value != None and isinstance(temp_out_dir, str):
                        command = ret_value.replace(placeholder, temp_out_dir)
            else:
                logger.error("Function %s does not return list or string object"
                             % func.__name__)
                sys.exit(1)
            return(command)
        return(inner)

    def get_run(self):
        if isinstance(self._eop, pipeline_info.PipelineInfo):
            run_info = self._eop.get_exec_group().get_run()
        elif isinstance(self._eop, exec_group.ExecGroup):
            run_info = self._eop.get_run()
        else:
            run_info = None
        return(run_info)

    def set_command(self, command):
        if not isinstance(command, list):
            logger.error("Given non-list command: %s" % command)
            sys.exit(1)
        self._command = command

    def set_stdout_path(self, stdout_path):
        if stdout_path != None:
            self._stdout_path = stdout_path

    @replace_output_dir_du_jour
    def get_stdout_path(self):
        return(self._stdout_path)

    def set_stderr_path(self, stderr_path):
        if stderr_path != None:
            self._stderr_path = stderr_path

    @replace_output_dir_du_jour
    def get_stderr_path(self):
        return(self._stderr_path)

    @replace_output_dir_du_jour
    def get_command(self):
        '''
        Return command after replacing all file inside the destination
        directory with relative paths.
        '''
        cmd = self._command
        run = self.get_run()
        working_dir = run.get_temp_output_directory()
        destination = run.get_step().get_pipeline().config['destination_path']
        diff = os.path.relpath(destination, working_dir)
        cmd = [arg.replace(destination, diff) for arg in cmd]
        return(cmd)
