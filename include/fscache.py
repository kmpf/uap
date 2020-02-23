import os
import yaml

class FSCache:
    '''
    Use this class if you expect to make the same os.path.* calls many
    times during a short time. The first time you call a method with certain
    arguments, the call is made, but all subsequent calls are served from a
    cache.

    Usage example::

        # Instantiate a new file system cache.
        fsc = FSCache()

        # This call will stat the file system.
        print(fsc.exists('/home'))

        # This call will leave the file system alone, the cached result will be returned.
        print(fsc.exists('/home'))

    You may call any method which is available in os.path.
    '''

    def __init__(self):
        self.cache = dict()

    def load_yaml_from_file(self, path):
        if not 'load_yaml_from_file' in self.cache:
            self.cache['load_yaml_from_file'] = dict()

        if path in self.cache['load_yaml_from_file']:
            return self.cache['load_yaml_from_file'][path]

        f = open(path, 'r')
        data = yaml.load(f, Loader=yaml.FullLoader)
        f.close()
        self.cache['load_yaml_from_file'][path] = data
        return data

    def clear(self):
        self.cache = dict()

    def __getattr__(self, name):

        def method(*args):

            # if the function was already called with the same args, return
            # the result from the cache
            if name in self.cache:
                if args in self.cache[name]:
                    return self.cache[name][args]

            # otherwise, make the call and store the result in the cache
            try:
                result = getattr(os.path, name)(*args)
                if not name in self.cache:
                    self.cache[name] = dict()
                self.cache[name][args] = result
                return result
            except AttributeError:
                raise AttributeError("Module os.path has no method '%s'" % name)

        return method
