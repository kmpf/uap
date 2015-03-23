import sys
import re
import shutil
import yaml
from yaml import dump 
from functools import wraps

from abstract_step import *
import command
import pipeline
import process_pool
import run as run_module

#class Step_info(object):
#    '''
#    Wrapper that holds all infos about how a step needs to be configured
#    '''
#    def __init__(self, abstract_step):
#        self._abstract_step = abstract_step
#        self._runs = dict()
#        
#    def new_run(self, run_name):
#        run = Run_info(self, run_name)
#        self._runs[run_name] = run
#        return run
#
#    def get_run(self, run_name):
#        return self._runs[run_name]
#
#    def get_runs(self):
#        for run in self._runs:
#            yield self._runs[run]
#
#    def get_abstract_step(self):
#        return self._abstract_step
#
#    def get_output_directory_du_jour_placeholder(self):
#        '''
#        Returns a placeholder for the temporary output directory, which
#        needs to be replaced by the actual temp directory inside the
#        execute() method
#        '''
#        return("<%s-output-directory-du-jour>" %  
#               str(self._abstract_step.__class__.__name__))
#
#class Run_info(object):
#    '''
#    Wrapper that holds all information about a single run
#    '''
#    def __init__(self, step, name):
#        self._step = step
#        self._name = str
#        self._exec_groups = list()
#        self._public_info = dict()
#        self._in_out_connection = list()
#        self._connections = list()
#        self._output_files = list()
#        self._input_files = list()
#        self._run_module = None
#        self.set_name(name)

#    def new_exec_group(self):
#        exec_group = Exec_group(self)
#        self._exec_groups.append(exec_group)
#        return exec_group
        
#    def get_exec_groups(self):
#        return self._exec_groups

#    def get_step(self):
#        return self._step

#    def set_name(self, name):
#        self._name = name
        
#    def get_name(self):
#        return self._name

#    def get_output_files_only(self):
#        return self._output_files
    
#    def set_run_module(self, run_object):
#        if isinstance(run_object, run_module.Run):
#            self._run_module = run_object
#            self._add_output_files_to_run_module()
#            self._add_public_info_to_run_module()
#            self._store_step_info_in_run_module()
#        else:
#            raise StandardError("Received object of type '%s' expected type "
#                                "'%s'" %
#                                (run_object.__class__.__name__, 
#                                 run_module.Run.__name__))

#    def get_run_module(self):
#        return self._run_module

#    def set_public_info(self, **kwargs):
#        for k, v in kwargs.items():
#            self._public_info[k] = v
            
#    def get_complete_public_info(self):
#        return self._public_info
        
#    def _add_public_info_to_run_module(self):
#        '''
#        Add public info to run_module.Run object
#        '''
#        if isinstance(self._run_module, run_module.Run):
#            # add public information if need be
#            for info_key, info_value in self.get_complete_public_info().items():
#                run.add_public_info(info_key, info_value)
#        else:
#            raise StandardError("Call set_run(run_object) prior to calling "
#                                "_add_public_info_to_run_module()")
        
#    def add_empty_output_connection(self, connection):
#        self._in_out_connection.append( (connection, None, None) )
#        
#    def add_output_file(self, connection, output_file, input_files):
#        self._connections.append(connection)
#        self._output_files.append(output_file)
#        self._input_files.append(input_files)
#        self._in_out_connection.append( (connection,
#            output_file, input_files) )
#        return "%s/%s" % (self._step.get_output_directory_du_jour_placeholder(),
#                          output_file)
#        
#    def _add_output_files_to_run_module(self):
#        '''
#        Add output files to run_module.Run object
#        '''
#        if isinstance(self._run_module, run_module.Run):
#            # ioc stands for in_out_connection
#            for ioc in self.get_output_file_triple():
#                if ioc[1] == None and ioc[2] == None:
#                    self._run_module.add_empty_output_connection(ioc[0])
#                else:
#                    self._run_module.add_output_file(ioc[0], ioc[1], ioc[2])
#        else:
#            raise StandardError("Call set_run(run_object) prior to calling "
#                                "_add_output_files_to_run_module()")
#
#    def _store_step_info_in_run_module(self):
#        '''
#        Store step_info in run_module.Run.add_private_info()
#        '''
#        if isinstance(self._run_module, run_module.Run):
#            self._run_module.add_private_info('step_info', self._step)
#
#    def get_output_file_triple(self):
#        return self._in_out_connection
#
#class Exec_group(object):
#    def __init__(self, run_info):
#        self._run_info = run_info
#        self._pipes_and_commands = list()
#        
#    def new_pipeline(self):
#        pipeline = Pipeline_info(self)
#        self._pipes_and_commands.append(pipeline)
#        return pipeline
#
#    def new_command(self, command, stderr_path=None, stdout_path=None):
#        command = Command_info(self, command, stderr_path, stdout_path)
#        self._pipes_and_commands.append(command)
#        return command
#
#    def get_pipes_and_commands(self):
#        return self._pipes_and_commands
#
#    def get_run_info(self):
#        return self._run_info
#
#class PipelineInfo(object):
#    def __init__(self, exec_group):
#        self._exec_group = exec_group
#        self._commands = list()
#
#    def new_command(self):
#        command = Command_info(self)
#        self._commands.append(command)
#        return command
#        
#    def get_commands(self):
#        return self._commands
#
#    def get_exec_group(self):
#        return self._exec_group
#
#class Command_info(object):
#    def __init__(self, eop, command, stdout_path=None, stderr_path=None):
#        # eop = exec_group or pipeline
#        self._eop = eop
#        self._command = list()
#        self._stdout_path = stdout_path
#        self._stderr_path = stderr_path
#        self._tool = str
#        self._output_files_per_connection = dict()
#
#        for _ in command:
#            if not isinstance(_, str):
#                raise StandardError("Non-string element %s in command %s" 
#                                    % (_, command))
#            self._command.append(_)
#        # get tool needed
#        self._tool = command[0]
#
#    def replace_output_dir_du_jour(func):
#        def inner(self, *args):
#            run_info = None
#            if isinstance(self._eop, Pipeline_info):
#                run_info = self._eop.get_exec_group().get_run_info()
#            elif isinstance(self._eop, Exec_group):
#                run_info = self._eop.get_run_info()
#            # Collect info to replace du_jour placeholder with temp_out_dir
#            step = run_info.get_step()
#            placeholder = step.get_output_directory_du_jour_placeholder()
#            temp_out_dir = step.get_abstract_step().get_output_directory_du_jour()
#
#            command = list()
#            ret_value = func(*args)
#            if isinstance(ret_value, list):
#                for string in list_of_strings:
#                    if string != None and placeholder in string:
#                        command.append(
#                            string.replace(placeholder, temp_out_dir))
#            elif isinstance(ret_value, str):
#                if ret_value != None and placeholder in ret_value:
#                        command.append(
#                            ret_value.replace(placeholder, temp_out_dir))
#            else:
#                raise StandardException("Function %s does not return list or "
#                                        "string object" % 
#                                        func.__class__.__name__)
#                
#            return(command)
#        return inner
#
#
#    def set_command(self, command):
#        if not isinstance(command, list):
#            raise StandardError("Given non-list command: %s" % command)
#        self._command = command
#
#    def set_stdout_path(self, stdout_path):
#        if stdout_path != None:
#            self._stdout_path = stdout_path
#
#    @replace_output_dir_du_jour
#    def get_stdout_path(self):
#        stdout_path = self._replace_output_dir_du_jour(self._stdout_path)
#        print(stdout_path)
#        return(stdout_path)
#
#    def set_stderr_path(self, stderr_path):
#        if stderr_path != None:
#            self._stderr_path = stderr_path
#
#    @replace_output_dir_du_jour
#    def get_stderr_path(self):
#        stderr_path = self._replace_output_dir_du_jour(self._stderr_path)
#        print(stderr_path)
#        return(stderr_path)
#
#    @replace_output_dir_du_jour
#    def get_command(self):
#
#        command = list()
#        for _ in self._command:
#            _ = self._replace_output_dir_du_jour(_)
#            command.append(_)
#        print(command)
#        
#        return(command)
#

        

class Fastqc(AbstractStep):
    '''
    | The fastqc step  is a wrapper for the fastqc tool. 
    | It generates some quality metrics for fastq files.
    | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
    | For this specific instance only the zip archive is preserved
    '''
    
    def __init__(self, pipeline):
        super(Fastqc, self).__init__(pipeline)
        
        self.set_cores(1) # muss auch in den Decorator
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read_fastqc_report')
        self.add_connection('out/first_read_log_stderr')
        self.add_connection('out/second_read_fastqc_report')
        self.add_connection('out/second_read_log_stderr')
        
        self.require_tool('fastqc') # auch in den Decorator?
        self.add_option('contaminent-file', str, optional =True) # muss seperat bleiben sind nicht nur Optionen fuer tools


    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for read in read_types:
            for run_id in run_ids_connections_files.keys():
                connection = 'in/%s' % read
                input_paths = run_ids_connections_files[run_id][connection]

                run = self.new_run(run_id)
                if input_paths == [None]:
                    run.add_empty_output_connection("%s_fastqc_report" % read)
                    run.add_empty_output_connection("%s_log_stderr" % read)
                else:
                    fastqc_exec_group = run.new_exec_group()
                    fastqc = ['fastqc',
                              '--noextract', '-o',
                              self.get_output_directory_du_jour_placeholder()]
                    fastqc.extend(input_paths)
                    
                    fastqc_command = fastqc_exec_group.new_command(
                        fastqc,
                        stderr_path = run.add_output_file(
                            "%s_log_stderr" % read, 
                            "%s%s-fastqc-log_stderr.txt" % 
                            (run_id, read_types[read]),
                            input_paths))
                    
                    mv_exec_group = run.new_exec_group()
                    input_base = os.path.basename(
                        input_paths[0]).split('.', 1)[0]
                    mv = ['mv',
                          os.path.join(self.get_output_directory_du_jour_placeholder(),
                                       ''.join([input_base, '_fastqc.zip'])),
                          run.add_output_file("%s_fastqc_report" % read,
                            "%s%s-fastqc.zip" % (run_id, read_types[read]),
                            input_paths)]
                    mv_command = mv_exec_group.new_command(mv)

    def declare_runs(self):
        # Was muss hier alles passieren damit es funktioniert?
        # * es muessen alle runs definiert werden
        # * pro run muessen alle public/private Infos gesetzt werden
        # * es MUESSEN die Output Dateien den Connections zugeordnet werden

        # fetch all incoming run IDs which produce reads...

        in_connections = self.get_in_connections()
        run_ids_connections_files = dict()
        for in_connection in in_connections:
            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection(in_connection):
                # das macht den schoenen Generator kaputt den Micha mal gebaut hat
                if not run_id in run_ids_connections_files:
                    run_ids_connections_files[run_id] = dict()
                if not in_connection in run_ids_connections_files[run_id]:
                    run_ids_connections_files[run_id][in_connection] = input_paths
        
        self.runs(run_ids_connections_files)

        # inspect each run
#        for run_info in self.get_runs():
#            # declare run here
#            with self.declare_run(run_info.get_run_id()) as run:
#                # make run object known to run_info
#                run_info.set_run_module(run)
                
### Finished ###

    def execute(self, run_id, run):
        # Ich muss noch ne Loesung finden um hier beliebigen Python Code auszufuehren
        # exec() oder eval()?
        # get the step_info object back
#        step = run.get_private_info('step')

        # get the temporary output directory info
#        temp_out_dir = self.get_output_directory_du_jour()
#        placeholder = step_info.get_output_directory_du_jour()
#        def fix_du_jour_issue(in_command):
#            out_command = list()
#            for _ in in_command:
#                out_command.append(_.replace(placeholder, temp_out_dir))
#            return out_command
#
        # get run_info objects
        run_info = self.get_run(run_id)
        
        # for each exec_group in that run ...
        for exec_group in run_info.get_exec_groups():
            # ... create a process pool
            with process_pool.ProcessPool(self) as pool:
                for poc in exec_group.get_pipes_and_commands():
                    # for each pipe or command (poc)
                    # check if it is a pipeline ...
                    if isinstance(poc, Pipeline_info):
                        # ... create a pipeline ...
                        with pool.Pipeline(pool) as pipeline:
                            for command in poc.get_commands:
#                                command.set_command(
#                                    fix_du_jour_issue(command.get_command())
#                                )
                                pipeline.append(
                                    command.get_command(),
                                    stdout_path = command.get_stdout_path(),
                                    stderr_path = command.get_stderr_path())
                    elif isinstance(poc, CommandInfo):
#                        poc.set_command( 
#                            fix_du_jour_issue(poc.get_command())
#                        )
                        print(poc.get_command())
                        pool.launch(
                            poc.get_command(),
                            stdout_path = poc.get_stdout_path(),
                            stderr_path = poc.get_stderr_path())
                    else:
                        raise StandardError("[%s]" % self.__class__.__name__)
