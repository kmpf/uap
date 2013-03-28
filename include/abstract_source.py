import sys
sys.path.append('./include/sources')
import inspect

def get_source_class_for_key(key):
    classes = [_ for _ in inspect.getmembers(__import__(key), inspect.isclass) if AbstractSource in _[1].__bases__]
    if len(classes) != 1:
        raise StandardError("need exactly one subclass of AbstractSource in " + key)
    return classes[0][1]

class AbstractSource(object):
    def __init__(self, pipeline):
        self.pipeline = pipeline
