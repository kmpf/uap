#!/usr/bin/env python

import ast
import inspect
import glob
import os
import string
import sys
from uaperrors import UAPError
from logging import getLogger
logger = getLogger('uap_logger')

from abstract_step import AbstractSourceStep, AbstractStep
from pipeline import Pipeline

def main(args):

    # Define Class to return all 'self.add_connection' calls
    class AddConnectionLister(ast.NodeVisitor):
        def visit_Call(self, node):
            value, attr = (str(), str())
            try:
                value = node.func.value.id
                attr = node.func.attr
            except AttributeError as e:
                pass
            if value == 'self' and attr == 'add_connection':
                for s in node.args:
                    print("  - '%s'" % s.s)
            self.generic_visit(node)


    # Define Class to return all 'self.require_tool' calls
    class RequireToolLister(ast.NodeVisitor):
        def visit_Call(self, node):
            value, attr = (str(), str())
            try:
                value = node.func.value.id
                attr = node.func.attr
            except AttributeError as e:
                pass
            if value == 'self' and attr == 'require_tool':
                for s in node.args:
                    print("  - '%s'" % s.s)

            self.generic_visit(node)
            
    # Define Class to return all 'self.add_option' calls
    class AddOptionLister(ast.NodeVisitor):
        def visit_Call(self, node):
            value, attr = (str(), str())
            try:
                value = node.func.value.id
                attr = node.func.attr
            except AttributeError as e:
                pass
            if value == 'self' and attr == 'add_option':
                try:
                    print("   - '%s' (expects %s)" % (node.args[0].s, node.args[1].id) )
                except:
                    pass

                for k in node.keywords:
                    try:
                        print("      |_ %s=%s " % (k.arg, k.value.id))
                        continue
                    except:
                        pass
                    try:
                        print("      |_ %s=%s " % (k.arg, k.value.s))
                        continue
                    except:
                        pass
                    try:
                        print("      |_ %s=%s " % (k.arg, ", ".join(
                            [x.s for x in k.value.elts])))
                    except:
                        pass
                        

            self.generic_visit(node)

    def is_key_a_step(key, step_type, state=False):
        '''
        Check if given key belongs to a loadable class
        '''
        #been there 
        if state == True:
            return False 

        res = False 
        for name, cl in inspect.getmembers(__import__(key), inspect.isclass):
                if  cl.__module__ == key:
                    if issubclass(cl, step_type):
                        res = True
        return res
        

    steps_path = os.path.dirname(os.path.realpath(__file__))

    # Assemble list of all files in ../sources and ../steps which end on .py
    # Check if they are loadable objects of type AbstractStep or AbstractSourceStep

    if args.step:
        is_step = False
        for cl in [AbstractSourceStep, AbstractStep]:
            if is_key_a_step(args.step, cl, is_step):
                is_step = True
                step_file = '../sources/%s.py' % args.step \
                        if cl == AbstractSourceStep else \
                           '../steps/%s.py' % args.step

                step_file = os.path.join(steps_path, step_file)
                
                if not os.path.exists(step_file):
                    raise UAPError("Step file %s does not exists." %
                                 step_file)
                
                with open(step_file) as f:
                    fc = f.read()
                    
                tree = ast.parse(fc)
                    
                step = AbstractStep.get_step_class_for_key(args.step)
                print("Step: %s" % args.step)
                print("* General Information:")
                print(step.__doc__)
                print("* Provided Connections:")
                AddConnectionLister().visit(tree)                
                print("* Required Tools:")
                RequireToolLister().visit(tree)
                print("* Available Options:")
                AddOptionLister().visit(tree)
                #print(ast.dump(tree) )
                #print(rt.body[0].value.args)
                
        if not is_step:
            raise UAPError("'%s' is neither a source nor a processing step"
                         % args.step)

    else:
        source_steps = sorted([
            os.path.splitext( os.path.basename(f) )[0] \
            for f in glob.glob(os.path.join(steps_path, '../sources/*.py'))
        ])

        proc_steps = sorted([
            os.path.splitext( os.path.basename(f) )[0] \
            for f in glob.glob(os.path.join(steps_path, '../steps/*.py'))
        ])

        print("\nAvailable steps (Ordered by type):")
        print("==================================")
        print("Source steps:")
        print("-------------")
        for s in source_steps:
            if is_key_a_step(s, AbstractSourceStep):
#            if is_key_a_step(s, AbstractStep):
                print("- %s" % s)

        print("\nProcessing steps:")
        print("-----------------")
        for s in proc_steps:
            if is_key_a_step(s, AbstractStep):
                print("- %s" % s)
