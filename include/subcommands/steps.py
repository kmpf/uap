#!/usr/bin/env python

import inspect
import glob
import os
import string
import sys
import yaml
from collections import OrderedDict
import re

from uaperrors import UAPError
from logging import getLogger
logger = getLogger('uap_logger')

from abstract_step import AbstractSourceStep, AbstractStep
from pipeline import Pipeline

def main(args):

    steps_path = os.path.dirname(os.path.realpath(__file__))

    # Assemble list of all files in ../sources and ../steps which end on .py
    # Check if they are loadable objects of type AbstractStep or AbstractSourceStep

    if args.step:
        is_step = False
        for cl in [AbstractSourceStep, AbstractStep]:
            if is_key_a_step(args.step, cl, is_step):
                is_step = True

                step_class = AbstractStep.get_step_class_for_key(args.step)
                step = step_class(None)
                print('General Information:')
                print(step.__doc__)
                for option in step._defined_options.values():
                    option['types'] = option['types string']
                    del option['types string']
                report = OrderedDict([
                    ('In Connections', connections_dict(step, 'in')),
                    ('Out Connections', connections_dict(step, 'out')),
                    ('Required Tools', step._tools.keys()),
                    ('Available Options', step._defined_options),
                    ('CPU Cores', step.get_cores())
                ])
                print(yaml.dump(report, Dumper=MyDumper, default_flow_style=False))

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
                if args.details:
                    print(s)
                    print('_'*len(s))
                    step_class = AbstractStep.get_step_class_for_key(s)
                    print(step_class.__doc__)
                else:
                    print("- %s" % s)

        print("\nProcessing steps:")
        print("-----------------")
        for s in proc_steps:
            if is_key_a_step(s, AbstractStep):
                if args.details:
                    print(s)
                    print('-'*len(s))
                    step_class = AbstractStep.get_step_class_for_key(s)
                    print(step_class.__doc__)
                else:
                    print("- %s" % s)

class literal(str):
    pass

def literal_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
yaml.add_representer(literal, literal_presenter)

def ordered_dict_presenter(dumper, data):
    return dumper.represent_dict(data.items())
yaml.add_representer(OrderedDict, ordered_dict_presenter)

class MyDumper(yaml.Dumper):

    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)

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


def connections_dict(step, way, strip_prefix=True):
    """
    Returns a dict of connection information for way = 'in' or 'out'.
    """
    if way == 'in':
        connections = step.get_in_connections(strip_prefix=strip_prefix)
    elif way == 'out':
        connections = step.get_out_connections(strip_prefix=strip_prefix)
    else:
        raise ValueError('The argument `way` needs to be "in" or "out".')
    report = dict()
    for conn in connections:
        report[conn] = dict()
        report[conn]['optional'] = conn in step._optional_connections
        if conn in step._connection_formats.keys():
            report[conn]['format'] = step._connection_formats[conn]
        wconn = way + '/' + conn
        if wconn in step._connection_descriptions.keys():
            report[conn]['description'] = step._connection_descriptions[wconn]
    return report
