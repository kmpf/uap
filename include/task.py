import copy
import csv
import datetime
import glob
import os
import re
import StringIO
import sys
import yaml

import abstract_step
import pipeline

class Task(object):

    def __init__(self, pipeline, step, run_id):
        self.pipeline = pipeline
        self.step = step
        self.run_id = run_id

    def __str__(self):
        task_state = self.step.get_run_state(self.run_id)
        return '[' + task_state[0].lower() + '] ' + str(self.step.get_step_id()) + '/' + self.run_id

    def get_task_state(self):
        return self.step.get_run_state(self.run_id)

    def run(self):
        task_state = self.step.get_run_state(self.run_id)
        if task_state != self.pipeline.states.READY:
            raise StandardError(str(self) + ' cannot be run yet.')
        self.step.run(self.run_id)
