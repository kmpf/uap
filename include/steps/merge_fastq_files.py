"""

"""


from abstract_step import *


class MergeFastqFiles(AbstractStep):
    """Merge all fastq files of a sample.

        
    """
    
    def __init__(self, pipeline):
        super(MergeFastqFiles, self).__init__(pipeline)
        
        self.set_cores(1) # muss auch in den Decorator
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        
        self.require_tool('cat')

    def runs(self, run_ids_connections_files):
        '''
        self.runs() should be a replacement for declare_runs() and execute_runs()
        All information given here should end up in the step object which is 
        provided to this method.
        '''
        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            run = self.new_run(run_id)
            for read in read_types:
                connection = 'in/%s' % read
                input_paths = run_ids_connections_files[run_id][connection]

                if input_paths == [None]:
                    run.add_empty_output_connection("%s" % read)
                else:
                    cat_exec_group = run.new_exec_group()
                    cat = ['cat']
                    cat.extend(input_paths)
                    
                    cat_command = cat_exec_group.add_command(
                        cat,
                        stdout_path = run.add_output_file(
                            "%s" % read, 
                            "%s%s.fastq" % 
                            (run_id, read_types[read]),
                            input_paths))
                    

