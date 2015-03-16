import sys
from abstract_step import *
import pipeline
import re
import process_pool
import shutil
import yaml
from yaml import dump 

from functools import wraps


class Step_info(object):
    '''
    Wrapper that holds all infos about how a step needs to be configured
    '''
    def __init__(self):
        self._runs = dict()
        
    def new_run(self, run_name):
        run = Run_info(run_name)
        self._runs[run_name] = run
        return run

    def get_run(self, run_name):
        return self._runs[run_name]

    def get_runs(self):
        for run in self._runs:
            yield self._runs[run]

    def get_output_directory_du_jour(self, step):
        '''
        Returns a placeholder for the temporary output directory, which
        needs to be replaced by the actual temp directory inside the
        execute() method
        '''
        return "<%s-output-directory-du-jour>" % str(step.__class__.__name__)

class Run_info(object):
    '''
    Wrapper that holds all information about a single run
    '''
    def __init__(self, name):
        self._name = str
        self._exec_groups = list()
        self._public_info = dict()
        self._in_out_connection = list()
        self.set_name(name)

    def new_exec_group(self):
        exec_group = Exec_group()
        self._exec_groups.append(exec_group)
        return exec_group
        
    def get_exec_groups(self):
        return self._exec_groups

    def set_name(self, name):
        self._name = name
        
    def get_name(self):
        return self._name

    def info(self, **kwargs):
        for k, v in kwargs.items():
            self._public_info[k] = v
            
    def get_info(self):
        return self._public_info
        
    def add_output_file(self, connection, output_file, input_files):
        self._in_out_connection.append( (connection, 
            output_file, input_files) )
        return output_file
            
    def get_output_file(self):
        return self._in_out_connection

    def add_empty_output_connection(self, connection):
        self._in_out_connection.append( (connection, None, None) )
        
class Exec_group(object):
    def __init__(self):
        self._pipes_and_commands = list()
        
    def new_pipeline(self):
        pipeline = Pipeline_info()
        self._pipes_and_commands.append(pipeline)
        return pipeline

    def new_command(self, command, stderr_path=None, stdout_path=None):
        command = Command_info(command, stderr_path, stdout_path)
        self._pipes_and_commands.append(command)
        return command

    def get_pipes_and_commands(self):
        return self._pipes_and_commands

class Pipeline_info(object):
    def __init__(self):
        self._commands = list()

    def new_command(self):
        command = Command_info()
        self._commands.append(command)
        return command
        
    def get_commands(self):
        return self._commands

class Command_info(object):
    def __init__(self, command, stdout_path=None, stderr_path=None):

        self._command = list()
        self._stdout_path = stdout_path
        self._stderr_path = stderr_path
        self._tool = str
        self._output_files_per_connection = dict()

        for _ in command:
            if not isinstance(_, str):
                raise StandardError("Non-string element %s in command %s" 
                                    % (_, command))
            print(_)
            self._command.append(_)
        # get tool needed
        self._tool = command[0]

    def get_command(self):
        return self._command

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


    def runs(self, run_ids_connections_files, step):
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

                run = step.new_run(run_id)
                if input_paths == [None]:
                    run.add_empty_output_connection("%s_fastqc_report" % read)
                    run.add_empty_output_connection("%s_log_stderr" % read)
                else:
                    fastqc_exec_group = run.new_exec_group()
                    fastqc = ['fastqc',
                              '--noextract', '-o',
                              step.get_output_directory_du_jour(self)]
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
                          os.path.join(step.get_output_directory_du_jour(self),
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
        
        step_info = Step_info()
        self.runs(run_ids_connections_files, step_info)

        # inspect each run
        for run_info in step_info.get_runs():
            # declare run here
            with self.declare_run(run_info.get_name()) as run:
                # add public information if need be
                for info_key, info_value in run_info.get_info().items():
                    run.add_public_info(info_key, info_value)
                    
                # setup all connections
                # ioc stands for in_out_connection
                for ioc in run_info.get_output_file():
                    if ioc[1] == None and ioc[2] == None:
                        run.add_empty_output_connection(ioc[0])
                    else:
                        run.add_output_file(ioc[0], ioc[1], ioc[2])
                    
                # Finally add step object to run
                run.add_private_info('step_info', step_info)
                
### Finished ###

    def execute(self, run_id, run):
        # Ich muss noch ne Loesung finden um hier beliebigen Python Code auszufuehren
        # exec() oder eval()?

        # get the temporary output directory info
        temp_out_dir = self.get_output_directory_du_jour()

        # get the step_info object back
        step_info = run.get_private_info('step_info')

        # iterate over all run_info objects
        for run_info in step_info.get_runs():

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
                                    print(command.get_command())
                                    pipeline.append(
                                        command.get_command(),
                                        stdout_path = command.get_stdout_path(),
                                        stderr_path = command.get_stderr_path())
                        elif isinstance(poc, Command_info):
                            print(poc.get_command())
                            pool.launch(
                                poc.get_command(),
                                stdout_path = command.get_stdout_path(),
                                stderr_path = command.get_stderr_path())
                        else:
                            raise StandardError("[%s]" % self.__class__.__name__)
