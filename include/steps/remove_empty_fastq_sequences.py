import sys
from abstract_step import *

class RemoveEmptyFastqSequences(AbstractStep):
    def __init__(self):
        super(RemoveEmptyFastqSequences, self).__init__()
