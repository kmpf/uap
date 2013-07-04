import abstract_step
import os
import yaml

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

    def get_task_state_basic(self):
        '''
        Proxy method for step.get_run_state().
        '''
        return self.step.get_run_state_basic(self.run_id)

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
            result.add(self.pipeline.task_for_task_id[self.pipeline.task_id_for_output_file[path]])
        return list(result)
    
    def volatilize_if_possible(self, srsly = False):
        result = set()
        if not self.step._options['_volatile']:
            return set()
        for path_a in self.pipeline.output_files_for_task_id[str(self)]:
            if abstract_step.AbstractStep.fsc.exists(path_a):
                # now check whether we can volatilize path A
                path_a_can_be_removed = True
                path_a_dependent_files = list()
                if path_a in self.pipeline.file_dependencies_reverse:
                    for path_b in self.pipeline.file_dependencies_reverse[path_a]:
                        path_a_dependent_files.append(path_b)
                        # don't check whether the output file B exists,
                        # it might also be volatile, rather check whether the
                        # task which creates B is finished
                        path_b_task = self.pipeline.task_for_task_id[self.pipeline.task_id_for_output_file[path_b]]
                        if path_b_task.get_task_state() != self.pipeline.states.FINISHED:
                            path_a_can_be_removed = False
                            break
                    
                if path_a_can_be_removed:
                    result.add(path_a)
                    if srsly:
                        print("Now volatilizing %s: %s" % (str(self), os.path.basename(path_a)))
                        info = dict()
                        info['self'] = dict()
                        info['self']['size'] = abstract_step.AbstractStep.fsc.getsize(path_a)
                        info['self']['mtime'] = abstract_step.AbstractStep.fsc.getmtime(path_a)
                        info['downstream'] = dict()
                        for path_b in path_a_dependent_files:
                            info['downstream'][path_b] = dict()
                            if abstract_step.AbstractStep.fsc.exists(path_b):
                                info['downstream'][path_b]['size'] = abstract_step.AbstractStep.fsc.getsize(path_b)
                                info['downstream'][path_b]['mtime'] = abstract_step.AbstractStep.fsc.getmtime(path_b)
                            else:
                                downstream_info = yaml.load(open(path_b + abstract_step.AbstractStep.VOLATILE_SUFFIX, 'r'))
                                info['downstream'][path_b]['size'] = downstream_info['self']['size']
                                info['downstream'][path_b]['mtime'] = downstream_info['self']['mtime']
                            
                        path_a_volatile = path_a + abstract_step.AbstractStep.VOLATILE_SUFFIX
                        with open(path_a_volatile, 'w') as f:
                            f.write(yaml.dump(info, default_flow_style = False))
                        
                        os.utime(path_a_volatile, (os.path.getatime(path_a), os.path.getmtime(path_a)))
                        os.unlink(path_a)
        return result
