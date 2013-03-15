import sys
sys.path.append('./include/steps')

def get_step_class_for_key(key):
    return __import__(key).Step()

class AbstractStep(object):
    def __init__(self):
        pass