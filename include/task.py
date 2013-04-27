class Task(object):
    '''
    A task represents a certain run of a certain step.
    '''

    def __init__(self, pipeline, step, run_id):
        self.pipeline = pipeline
        self.step = step
        self.run_id = run_id

    def __str__(self):
        return '%s/%s' % (self.step.get_step_id(), self.run_id)

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
        if task_state != self.pipeline.states.READY:
            raise StandardError("%s cannot be run yet." % self)
        self.step.run(self.run_id)
