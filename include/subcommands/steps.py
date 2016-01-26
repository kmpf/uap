#!/usr/bin/env python

import inspect
import glob
import os
import string

from abstract_step import AbstractSourceStep, AbstractStep

def main(args):

    def is_key_a_class(key):
        '''
        Check if given key belongs to a loadable class
        '''

        classes = [_ for _ in inspect.getmembers(__import__(key), 
                                                 inspect.isclass) \
                   if c in _[1].__bases__]


        for index, c in enumerate(check_classes):
            classes = [_ for _ in inspect.getmembers(__import__(key), 
                                                     inspect.isclass) \
                       if c in _[1].__bases__]
            print(classes[0][1])
            for k in range(index):
                classes = [_ for _ in classes if _[1] != check_classes[k]]
            if len(classes) > 0:
                if len(classes) != 1:
                    raise StandardError(
                        "need exactly one subclass of %s in %s" % (c, key))
                return classes[0][1]

        raise StandardError("No suitable class found for module %s." % key)


    # Assemble list of all files in ../sources and ../steps which end on .py
    # Check if they are loadable objects of type AbstractStep or AbstractSourceStep

    source_steps = [os.path.splitext( os.path.basename(f) )[0] \
                    for f in glob.glob('../sources/*.py')
                ]
    
    proc_steps = [os.path.splitext( os.path.basename(f) )[0] \
                    for f in glob.glob('../steps/*.py')
                ]

    print("Source steps:\n")
    for s in source_steps:
        is_key_a_class(s, AbstractSourceStep)

    print("Processing steps:\n")
    for s in proc_steps:
        is_key_a_class(s, AbstractStep)
