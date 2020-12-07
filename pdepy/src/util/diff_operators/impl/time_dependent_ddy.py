import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import time_dependent_op # diff_operators' core
import numpy as np
class td_ddy(time_dependent_op.time_dependent_operator):
    def __init__(self, dy, coefficient = lambda x, y: 1):
        implicit_stencil = None
        explicit_stencil = [((0, -1, -1), -1/(2*dy)), ((0, 1, -1), 1/(2*dy))]
        super().__init__(implicit_stencil, explicit_stencil=explicit_stencil, coefficient=coefficient)