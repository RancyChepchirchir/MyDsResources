import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import diff_op # diff_operators' core
import numpy as np
class ddy(diff_op.diff_operator):
    def __init__(self, dy, coefficient = lambda x, y: 1):
        super().__init__([((0, -1), -1/(2*dy)), ((0, 1), 1/(2*dy))], coefficient)