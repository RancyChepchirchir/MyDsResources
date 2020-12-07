import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import diff_op # diff_operators' core
import numpy as np
class laplacian2d(diff_op.diff_operator):
    def __init__(self, dx, dy, coefficient = lambda x, y: 1):
        super().__init__([((0, -1), 1/dy**2), ((-1, 0), 1/dx**2), ((0, 0), -2/dx**2-2/dy**2), ((1, 0), 1/dx**2), ((0, 1), 1/dy**2)], coefficient)