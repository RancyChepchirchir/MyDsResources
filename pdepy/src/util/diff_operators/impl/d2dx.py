import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import diff_op # diff_operators' core
from core.diff_op import irregular_stencil_node
import numpy as np
class d2dx(diff_op.diff_operator):
    def __init__(self, dx, coefficient = lambda x, y: 1):
        super().__init__([((-1, 0), 1/dx**2), ((0, 0), -2/dx**2), ((1, 0), 1/dx**2)], coefficient)

    def getIrregularStencil(self, dx, tau, direction):
        # direction = 1 represents positive y axis; -1 represents negative y axis
        c1, c2, c3, c4 = (tau-1)/(tau+2), 2*(2-tau)/(tau+1), (tau-3)/tau, 6/(tau*(tau+1)*(tau+2))
        if direction == diff_op.Direction.POSITIVE:
            return [irregular_stencil_node((-2, 0), c1/dx**2), irregular_stencil_node((-1, 0), c2/dx**2), irregular_stencil_node((0, 0), c3/dx**2), irregular_stencil_node((tau, 0), c4/dx**2, True)]
        elif direction == diff_op.Direction.NEGATIVE:
            return [irregular_stencil_node((-tau, 0), c4/dx**2, True), irregular_stencil_node((0, 0), c3/dx**2), irregular_stencil_node((1, 0), c2/dx**2), irregular_stencil_node((2, 0), c1/dx**2)]