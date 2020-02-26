from uaperrors import UAPError
import sys
import os
import yaml
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

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

    def get_pipeline(self):
        '''
        Returns the pipeline this task belongs to.
        '''
        return self.pipeline

    def get_step(self):
        '''
        Returns the step of this task.
        '''
        return self.step

    def get_run(self):
        '''
        Returns the run object for this task.
        '''
        return self.step.get_run(self.run_id)

    def get_task_state(self, do_hash=False):
        '''
        Proxy method for run.get_state().
        '''
        return self.get_run().get_state(do_hash=do_hash)

    def run(self):
        '''
        Run the task. Skip if it's already finished. Raise Exception
        if it's not ready.
        '''
        task_state = self.get_task_state()
        if task_state == self.pipeline.states.FINISHED:
            self.pipeline.notify("Skipping task: %s is already finished." % self)
            return
        if task_state == self.pipeline.states.WAITING:
            raise UAPError("%s cannot be run yet." % self)
        self.step.run(self.run_id)

    def generate_report(self):
        '''
        Generate the report for the task. Skip this if task is not finished yet.
        '''
        task_state = self.get_task_state()
        if task_state != self.pipeline.states.FINISHED:
            self.pipeline.notify("Skipping task: %s its not finished yet." % self)
            return
        self.step.generate_report(self.run_id)


    def input_files(self):
        '''
        Return a list of input files required by this task.
        '''
        result = set()
        run_info = self.get_run()
        for annotation, outfiles in run_info.get_output_files_abspath().items():
            for outpath, infiles in outfiles.items():
                if infiles != None:
                    for path in infiles:
                        result.add(path)
        return list(result)

    def output_files(self):
        '''
        Return a list of output files produced by this task.
        '''
        result = []
        run_info = self.get_run()
        for annotation, outfiles in run_info.get_output_files_abspath().items():
            for path in outfiles.keys():
                result.append(path)
        return result

    def get_parent_tasks(self):
        '''
        Returns a list of parent tasks which this task depends on.
        '''
        result = set()
        for path in self.input_files():
            result.add(self.pipeline.get_task_for_file(path))

        return list(result)

    def move_ping_file(self, bad_copy=True):
        '''
        Removes the queued ping of a task to remove the QUEUED state
        and optionally keeps "bad_copy" to mark the task as BAD.
        '''
        ping_file = self.get_run().get_queued_ping_file()
        self.step.move_ping_file(ping_file, bad_copy=bad_copy)

    def volatilize_if_possible(self, srsly = False):
        result = set()
        if not self.step.is_volatile():
            return result
        fsc = self.get_run().fsc
        for path_a in self.output_files():
            if not path_a: continue
            if fsc.exists(path_a):
                # now check whether we can volatilize path A
                path_a_can_be_removed = True
                path_a_dependent_files = list()
                if path_a in self.pipeline.file_dependencies_reverse:
                    for path_b in self.pipeline.file_dependencies_reverse[path_a]:
                        path_a_dependent_files.append(path_b)
                        # don't check whether the output file B exists,
                        # it might also be volatile, rather check whether the
                        # task which creates B is finished
                        path_b_task = self.pipeline.get_task_for_file(path_b)
                        if path_b_task.get_task_state() != \
                           self.pipeline.states.FINISHED:
                            path_a_can_be_removed = False
                            break

                if path_a_can_be_removed:
                    result.add(path_a)
                    if srsly:
                        self.pipeline.notify("Now volatilizing %s: %s" %
                              (str(self), os.path.basename(path_a)))
                        info = dict()
                        info['self'] = dict()
                        info['self']['size'] = fsc.getsize(path_a)
                        info['self']['mtime'] = fsc.getmtime(path_a)
                        info['downstream'] = dict()
                        for path_b in path_a_dependent_files:
                            info['downstream'][path_b] = dict()
                            if fsc.exists(path_b):
                                info['downstream'][path_b]['size'] = fsc.getsize(path_b)
                                info['downstream'][path_b]['mtime'] = fsc.getmtime(path_b)
                            else:
                                downstream_info = yaml.load(open(path_b + AbstractStep.VOLATILE_SUFFIX, 'r'), Loader=yaml.FullLoader)
                                info['downstream'][path_b]['size'] = downstream_info['self']['size']
                                info['downstream'][path_b]['mtime'] = downstream_info['self']['mtime']

                        path_a_volatile = path_a + AbstractStep.VOLATILE_SUFFIX
                        with open(path_a_volatile, 'w') as f:
                            f.write(yaml.dump(info, default_flow_style = False))

                        os.utime(path_a_volatile, (os.path.getatime(path_a), os.path.getmtime(path_a)))
                        os.unlink(path_a)
        return result
