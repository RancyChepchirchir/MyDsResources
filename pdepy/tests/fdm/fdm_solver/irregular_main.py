import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/fdm/solver/")
import fdm_solver as fdm
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/diff_operator_expression")
import diff_op_expression as expr
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/")
import diff_operators.impl.d2dx as d2dx, diff_operators.impl.d2dy as d2dy
import diff_operators.core.diff_op as diff_op
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/core/domain")
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/impl/dirichlet")
import dirichlet_bc as dr
import domain as dm
import unittest
import numpy as np
import math

def get_y_by_x(x, y):
    y1 = math.sqrt(1-x**2)
    y2 = -y1
    if abs(y-y2) < abs(y-y1):
        return (x, y2)
    else:
        return (x, y1)
def get_x_by_y(x, y):
    x1 = math.sqrt(1-y**2)
    x2 = -x1
    if abs(x-x2) < abs(x-x1):
        return (x2, y)
    else:
        return (x1, y)

class Test(unittest.TestCase):
    def setUp(self):
        domain = dm.domain(np.array([-1, -1]), np.array([1, 1]))
        inDomain = lambda x, y: x**2 + y**2 < 1
        onBoundary = lambda x, y: abs(x**2 + y**2 - 1) < np.spacing(1)
        getBoundaryValue = lambda x, y: 1
        getNearestPoint = [get_x_by_y, get_y_by_x]
        diri = dr.dirichlet_bc(inDomain, onBoundary, getBoundaryValue, domain, getNearestPoint)
        self.solver = fdm.fdm_solver(None, lambda x, y: 16*(x**2+y**2), diri)
    
    def test_solve_irregular(self):
        error = []
        for n in [9, 19, 39, 79]:
            dx, dy = 1/(n+1), 1/(n+1)
            a = d2dx.d2dx(dx)
            b = d2dy.d2dy(dy)
            expression = expr.diff_operator_expression([a, b])
            self.solver.diff_op_expression = expression
            u = self.solver.solve(2*n+1, 2*n+1)
            real_function = lambda x, y: (x**2+y**2)**2
            real_u = np.zeros(self.solver.vector_len)
            for i in range(self.solver.vector_len):
                real_u[i] = real_function(*self.solver._get_coord_by_offset(dx, dy, *self.solver.index_to_grid[i]))
            error.append(max(abs(u-real_u)))
        assert error[0] < 0.00996
        assert error[1] < 0.002497
        assert error[2] < 0.0006248
        assert error[3] < 0.000157

if __name__ == '__main__':
    unittest.main()