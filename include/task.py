class Task(object):
    '''
    A task represents a certain run of a certain step.
    '''

    def __init__(self, pipeline, step, run_id, run_index):
        self.pipeline = pipeline
        self.step = step
        self.run_id = run_id
        self.run_index = run_index

    def __str__(self):
        return '%s/%s' % (self.step.get_step_name(), self.run_id)

    def get_task_state(self):
        '''
        Proxy method for step.get_run_state().
        '''
        return self.step.get_run_state(self.run_id)

    def run(self):
        '''
        Run the task. Skip if it's already finished, otherwise raise 
        StandardError if it's not ready.
        '''
        task_state = self.get_task_state()
        if task_state == self.pipeline.states.FINISHED:
            print("Skipping task: %s is already finished." % self)
            return
        if task_state == self.pipeline.states.WAITING:
            raise StandardError("%s cannot be run yet." % self)
        self.step.run(self.run_id)

    def input_files(self):
        '''
        Return a list of input files required by this task.
        '''
        result = set()
        run_info = self.step.get_run_info()[self.run_id]
        for annotation, outfiles in run_info['output_files'].items():
            for outpath, infiles in outfiles.items():
                for path in infiles:
                    result.add(path)
        return list(result)

    def output_files(self):
        '''
        Return a list of output files produced by this task.
        '''
        result = []
        run_info = self.step.get_run_info()[self.run_id]
        for annotation, outfiles in run_info['output_files'].items():
            for path in outfiles.keys():
                result.append(path)
        return result

    def get_parent_tasks(self):
        '''
        Returns a list of parent tasks which this task depends on.
        '''
        result = set()
        for path in self.input_files():
            result.add(self.pipeline.task_for_output_file[path])
        return list(result)