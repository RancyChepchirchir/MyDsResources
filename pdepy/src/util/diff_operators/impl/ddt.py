import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import time_dependent_op # diff_operators' core
import numpy as np
class ddt(time_dependent_op.time_dependent_operator):
    def __init__(self, dt, coefficient = lambda x, y: 1):
        explicit_stencil = [((0, 0, -1), -1/dt), ((0, 0, 0), 1/dt)]
        super().__init__([((0, 0), 1/dt), ((0, -1), -1/dt)], explicit_stencil=explicit_stencil, coefficient=coefficient)