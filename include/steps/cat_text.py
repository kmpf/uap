import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class CatText(AbstractStep):
    '''
    cats text files together
    '''

    def __init__(self, pipeline):
        super(CatText, self).__init__(pipeline)
        
        self.set_cores(8)
        
        self.add_connection('in/text')
        self.add_connection('out/text')
        self.add_option('filenameEnding', str, optional=False)
        self.add_option('additionalFiles', list, optional=True)
        self.add_option('run_id', str, default='merged', optional=True)

        self.require_tool('cat')
#        self.require_tool('grep')
#        self.require_tool('sort')
#        self.require_tool('uniq')

    def runs(self, run_ids_connections_files):
        input_paths = []

        for run_id in run_ids_connections_files.keys():
            input_paths.extend(run_ids_connections_files[run_id]["in/text"])

        run_id = self.get_option('run_id')
        with self.declare_run(run_id) as run:

            if input_paths == [None]:
                run.add_empty_output_connection("text")
            

            out = run.add_output_file("text",
                                      "%s-combined.%s" %  
                                      (run_id, self.get_option('filenameEnding')),
                                       input_paths) 




#            out = run.add_output_file("text","foo", input_paths) 
                    
            with run.new_exec_group() as exec_group:
                cat = [self.get_tool('cat')]
                if self.is_option_set_in_config('additionalFiles'):
                    cat.extend(self.get_option('additionalFiles'))
                cat.extend(input_paths)
                exec_group.add_command(cat, stdout_path=out)




                            



