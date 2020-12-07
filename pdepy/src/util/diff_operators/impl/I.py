import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import diff_op # diff_operators' core
import numpy as np
class I(diff_op.diff_operator):
    def __init__(self, coefficient = lambda x, y: 1):
        super().__init__([((0, 0), 1)], coefficient)

    def get_implicit_stencil(self):
        return self.stencil
    
    def get_explicit_stencil(self):
        return self.stencil