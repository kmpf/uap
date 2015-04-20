import sys
from abstract_step import *
import pipeline
import re
import process_pool
import shutil
import yaml
from yaml import dump 

from functools import wraps

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
        
        self.add_connection('in/reads')
        self.add_connection('out/fastqc_report')
        self.add_connection('out/log_stderr')
        
        self.require_tool('fastqc') # auch in den Decorator
        self.add_option('contaminent-file', str, optional =True) # muss seperat bleiben sind nicht nur Optionen fuer tools

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

        def get_runs(self):
            return self._runs

    class Run_info(object):
        '''
        Wrapper that holds all information about a single run
        '''
        def __init__(self, name):
            self._name = str
            self._exec_groups = list()
            self._public_info = dict()
            self._connection_output_input = list()
            self.name(name)

        def new_exec_group(self):
            exec_group = Exec_group()
            self._exec_groups.append(exec_group)
            return exec_group
            
        def get_exec_groups(self):
            return self._exec_groups

        def name(self, name):
            self._name = name
            
        def get_name(self):
            return self._name

        def info(self, **kwargs):
            for k, v in kwargs.items():
                self._public_info[k] = v
                
        def get_info(self):
            return self._public_info
            
        def output_file(self, connection, output_file, input_files):
            self._connection_output_input.append( (connection, 
                output_file, input_files) )
                
        def get_output_file(self):
            return self._connection_output_input
            
    class Exec_group(object):
        def __init__(self):
            self._pipes_and_commands = list()
            
        def new_pipeline(self):
            pipeline = Pipeline_info()
            self._pipes_and_commands.append(pipeline)
            return pipeline

        def new_command(self):
            command = Command_info()
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
        def __init__(self, *args):

            self._command = list()
            self._stdout_path = None
            self._stderr_path = None
            self._tool = str
            self._output_files_per_connection = dict()

            for _ in args:
                if not isinstance(str, _):
                    raise StandardError("Non-string element %s in command %s" 
                                        % (_, args))
                self._command.append(_)
            # get tool needed
            self._tool = args[0]

        def get_command(self):
            return self._command


    def runs(self, run_ids_input_files, step):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        for run_id, input_paths in run_ids_input_files.items():
            # Combine or split input_paths here as necessary
            is_paired_end = self.find_upstream_info_for_input_paths(input_paths,
                'paired_end')
            first_read = self.find_upstream_info_for_input_paths(input_paths,
                'first_read')
            second_read = self.find_upstream_info_for_input_paths(input_paths,
                'second_read')

            # decide which read type we'll handle based on whether this is
            # paired end or not
            read_types = [first_read]
            if  is_paired_end:
                read_types.append(second_read)

            # put input files into R1/R2 bins (or one single R1 bin)
            input_path_bins = dict()
            for _ in read_types:
                input_path_bins[_] = list()

            for path in input_paths:
                which = misc.assign_string(os.path.basename(path), read_types)
                input_path_bins[which].append(path)

            # Add runs to step
            for which in read_types:
                run = step.new_run("%s%s" % (run_id, which))
                # set output files
                run.output_file("fastqc_report", "%s%s-fastqc.zip" 
                    % (run_id, which), input_path_bins[which])
                run.output_file("log_stderr", "%s%s-fastqc-log_stderr.txt" 
                    % (run_id, which), input_path_bins[which])

            # Hier definieren wir die ganzen Befehle die auszufuehren sind
            self.commands.append(
                [self.tool('fastqc'),
                 '--noextract', '-o',
                 self.out(fastqc_out_dir, )])
                      
        
    def declare_runs(self):
        # Was muss hier alles passieren damit es funktioniert?
        # * es muessen alle runs definiert werden
        # * pro run muessen alle public/private Infos gesetzt werden
        # * es MUESSEN die Output Dateien den Connections zugeordnet werden

        # fetch all incoming run IDs which produce reads...
        run_ids_input_files = dict()
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/reads'):
            # das macht den schoenen Generator kaputt den Micha mal gebaut hat
            run_ids_input_files[run_id] = input_paths
            
        step_info = Step_info()
        self.runs(run_ids_input_files, step_info)

        # inspect each run
        for run_info in step_info.get_runs():
            # declare run here
            with self.declare_run(run_info.get_name()) as run:
                # add public information if need be
                for info_key, info_value in run_info.get_info().items():
                    run.add_public_info(info_key, info_value)
                    
                # setup all connections
                for coi in run_info.get_output_file():
                # coi stands for connection_output_input
                    run.add_output_file(coi[0], coi[1], coi[2])

                # Finally add step object to run
                run.add_private_info('step_info', step_info)
                
### Finished ###

    def execute(self, run_id, run):
        # Ich muss noch ne Loesung finden um hier beliebigen Python Code auszufuehren
        # exec() oder eval()?

        # get the step_info object back
        step_info = run.get_private_info('step_info')
        # get the correct run_info object
        run_info = step_info.get_runs()[run_id]

        # for each exec_group in that run ...
        for exec_group in run_info.get_exec_groups():
            # ... create a process pool
            with process_pool.ProcessPool(self) as pool:
                # get the list of pipes and commands (pac)
                for pac in exec_group.get_pipes_and_commands():
                    # for each pipe or command (poc)
                    for poc in pac:
                        # check if it is a pipeline ...
                        if isinstance(Pipeline_info, poc):
                            # ... create a pipeline ...
                            with pool.Pipeline(pool) as pipeline:
                                for command in poc.get_commands:
                                    pipeline.append(command.get_command(),
                                        stdout_path = command.get_stdout_path(),
                                        stderr_path = command.get_stderr_path())


                # Fastqc only allows to write to a directory 
                fastqc_out_dir =  self.get_output_directory_du_jour()

                out_path = run.get_single_output_file_for_annotation('fastqc_report')
                in_path  = run.get_input_files_for_output_file(out_path)

                
                
                # set up processes                              
                fastqc = [self.get_tool('fastqc'), '--noextract', '-o', fastqc_out_dir]
                fastqc.extend(in_path)

                
                # create the pipeline and run it
                log_stderr = run.get_single_output_file_for_annotation('log_stderr')
                pipeline.append(fastqc, stderr_path = log_stderr)                   




        fastqc_default_name = run.get_private_info('fastqc_default_name')
        fastqc_report_basename  = fastqc_default_name + '.zip'

        full_path_zipped_fastqc_report = os.path.join(fastqc_out_dir,  fastqc_report_basename)


        
        try:
            os.rename(full_path_zipped_fastqc_report, out_path)
        except OSError:
            raise StandardError("os.rename failed of %s to %s" % full_path_zipped_fastqc_report, out_path) 


 
        # in case of:
        #unzipped_fastqc_report = os.path.join(fastqc_out_dir,  fastqc_default_name)
        #try:
        #    shutil.rmtree(unzipped_fastqc_report)
        #except OSError:
        #    raise StandardError('removing unzipped dir failes')
