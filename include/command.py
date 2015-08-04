"""

"""
import pipeline_info
import exec_group

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
            if not isinstance(_, str):
                raise StandardError("Non-string element %s in command %s" 
                                    % (_, command))
            self._command.append(_)
        # get tool needed
        self._tool = command[0]

    def replace_output_dir_du_jour(func):
        def inner(self, *args):
            run_info = None
            if isinstance(self._eop, pipeline_info.PipelineInfo):
                run_info = self._eop.get_exec_group().get_run()
            elif isinstance(self._eop, exec_group.ExecGroup):
                run_info = self._eop.get_run()
            # Collect info to replace du_jour placeholder with temp_out_dir
            step = run_info.get_step()
            placeholder = step.get_output_directory_du_jour_placeholder()
            temp_out_dir = step.get_output_directory_du_jour()
            
            command = None
            ret_value = func(self, *args)
            if isinstance(ret_value, list):
                command = list()
                for string in ret_value:
                    if string != None and placeholder in string:
                        command.append(
                            string.replace(placeholder, temp_out_dir))
                    else:
                        command.append(string)
            elif isinstance(ret_value, str):
                if ret_value != None and placeholder in ret_value:
                        command = ret_value.replace(placeholder, temp_out_dir)
            elif ret_value == None:
                command = None
            else:
                raise StandardError("Function %s does not return list or "
                                        "string object" % 
                                        func.__class__.__name__)
            return(command)
        return(inner)

    def set_command(self, command):
        if not isinstance(command, list):
            raise StandardError("Given non-list command: %s" % command)
        self._command = command

    def set_stdout_path(self, stdout_path):
        if stdout_path != None:
            self._stdout_path = stdout_path

    @replace_output_dir_du_jour
    def get_stdout_path(self):
#        stdout_path = self._replace_output_dir_du_jour(self._stdout_path)
        return(self._stdout_path)

    def set_stderr_path(self, stderr_path):
        if stderr_path != None:
            self._stderr_path = stderr_path

    @replace_output_dir_du_jour
    def get_stderr_path(self):
#        stderr_path = self._replace_output_dir_du_jour(self._stderr_path)
        return(self._stderr_path)

    @replace_output_dir_du_jour
    def get_command(self):

#        command = list()
#        for _ in self._command:
#            _ = self._replace_output_dir_du_jour(_)
#            command.append(_)
        return(self._command)
