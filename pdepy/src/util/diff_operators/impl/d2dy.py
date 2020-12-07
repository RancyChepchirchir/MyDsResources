import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/..')
from core import diff_op # diff_operators' core
from core.diff_op import irregular_stencil_node
import numpy as np
class d2dy(diff_op.diff_operator):
    def __init__(self, dy, coefficient = lambda x, y: 1):
        super().__init__([((0, -1), 1/dy**2), ((0, 0), -2/dy**2), ((0, 1), 1/dy**2)], coefficient)

    def getIrregularStencil(self, dy, tau, direction):
        # direction = 1 represents positive y axis; -1 represents negative y axis
        c1, c2, c3, c4 = (tau-1)/(tau+2), 2*(2-tau)/(tau+1), (tau-3)/tau, 6/(tau*(tau+1)*(tau+2))
        if direction == diff_op.Direction.POSITIVE:
            return [irregular_stencil_node((0, -2), c1/dy**2), irregular_stencil_node((0, -1), c2/dy**2), irregular_stencil_node((0, 0), c3/dy**2), irregular_stencil_node((0, tau), c4/dy**2, True)]
        elif direction == diff_op.Direction.NEGATIVE:
            return [irregular_stencil_node((0, -tau), c4/dy**2, True), irregular_stencil_node((0, 0), c3/dy**2), irregular_stencil_node((0, 1), c2/dy**2), irregular_stencil_node((0, 2), c1/dy**2)]