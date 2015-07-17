import sys
import abstract_step 
import copy
import misc
import os
import pipeline
import process_pool


class IOStep(abstract_step.AbstractStep):
    
    def __init__(self, pipeline, tool):
        super(IOStep, self).__init__(pipeline)
        
        self._tool = tool
        
        self.set_cores(4)
        self.add_connection('in/*')
        # It's not necessary to declare the "out/*" connection here, as the real
        # out connections are declare in declare_runs(). However, this is done
        # here for cosmetic purposes so that it shows up the step documentation.
        self.add_connection('out/*')
        
        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool(self._tool)
        

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/*'):
            with self.declare_run(run_id) as run:
                for in_path in input_paths:
                    annotation = self.get_annotation_for_input_file(in_path)
                    self.add_connection('out/%s' % annotation)
                    run.add_output_file(annotation, misc.append_suffix_to_path(os.path.basename(in_path), self._tool), [in_path])
                    run.new_exec_group()

    
    def execute(self, run_id, run):
        # process one file at a time
        for tag, output_file_info in run.get_output_files_abspath().items():
            for output_path, input_paths in output_file_info.items():
                with process_pool.ProcessPool(self) as pool:
                    with pool.Pipeline(pool) as pipeline:
                        cat4m = [self.get_tool('cat4m'), input_paths[0]]
                        pigz1 = [self.get_tool('pigz'), '--decompress', '--processes', '1', '--stdout']
                        tool_command = [self.get_tool(self._tool)]
                        tool_command.extend(self.tool_command_line())
                        pigz2 = [self.get_tool('pigz'), '--processes', '2', '--blocksize', '4096', '--stdout']
                
                        if output_path[-3:] == '.gz':
                            pipeline.append(cat4m)
                            pipeline.append(pigz1)
                            pipeline.append(tool_command)
                            pipeline.append(pigz2, stdout_path = output_path)
                        else:
                            pipeline.append(cat4m)
                            pipeline.append(tool_command, stdout_path = output_path)

    def tool_command_line(self):
        '''
        Raise NotImplementedError because every subclass must override this method.
        '''
        raise NotImplementedError()
