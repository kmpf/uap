import sys
sys.path.append('./include/steps')
import inspect

def get_step_class_for_key(key):
    classes = inspect.getmembers(__import__(key), inspect.isclass)
    #for c in classes:
        #print(c[1].__bases__)
        #if AbstractStep in c[1].__bases__:
            #print(c[0] + "ok")
        #print(inspect.getmro(c[1]))
    classes = [_ for _ in classes if AbstractStep in _[1].__bases__]
    if len(classes) != 1:
        raise StandardError("need exactly one subclass of AbstractStep in " + key)
    return classes[0][1]

class AbstractStep(object):
    def __init__(self):
        self.dependencies = []
        self.options = {}

    def set_options(self, options):
        self.options = options

    def add_dependency(self, parent):
        if not isinstance(parent, AbstractStep):
            raise StandardError("parent argument must be an AbstractStep")
        if parent == self:
            raise StandardError("cannot add a node as its own dependency")
        self.dependencies.append(parent)

    def __str__(self):
        s = self.__class__.__name__
        if len(self.options.keys()) > 0:
            s += " with options: " + str(self.options)
        if len(self.dependencies) > 0:
            s += " with dependencies: " + ', '.join(['(' + str(_) + ')' for _ in self.dependencies])
        return s

