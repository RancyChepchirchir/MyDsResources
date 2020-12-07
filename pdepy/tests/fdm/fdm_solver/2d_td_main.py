import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/fdm/solver/")
import fdm_solver as fdm
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/diff_operator_expression")
import diff_op_expression as expr
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/")
import diff_operators.impl.ddx as ddx, diff_operators.impl.ddy as ddy, diff_operators.impl.laplacian2d as laplacian2d,\
    diff_operators.impl.ddt as ddt, diff_operators.impl.time_dependent_d2dx as td_d2dx, diff_operators.impl.time_dependent_d2dy as td_d2dy
import diff_operators.core.diff_op as diff_op
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/core/domain")
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/impl/dirichlet")
import dirichlet_rectangle as dr
import domain as dm
import unittest
import numpy as np
import math


class Test(unittest.TestCase):
    def setUp(self):
        td_domain = dm.domain(np.array([0, 0]), np.array([1, 1]))
        td_inDomain = lambda x, y: 0 < x < 1 and 0 < y < 1
        td_onBoundary = lambda x, y: abs(x-1) < np.spacing(1) or abs(x) < np.spacing(1) \
            or abs(y) < np.spacing(1) or abs(y-1) < np.spacing(1)
        def td_getBV(x, y):
            if abs(x) < np.spacing(1) or abs(x-1) < np.spacing(1) or abs(y) < np.spacing(1) or abs(y-1) < np.spacing(1):
                return 0
            elif y <= 1/2 and (y <= x and y <= -x+1):
                return 20*y
            elif x <= 1/2 and (y >= x and y <= -x+1):
                return 20*x
            elif x >= 1/2 and (y <= x and y >= -x+1):
                return 20-20*x
            else:
                return 20-20*y
        f_s = lambda x, y: 0
        dirichlet = dr.dirichlet_rectangular_bc(td_inDomain, td_onBoundary, td_getBV, td_domain, 2)
        self.solver = fdm.fdm_solver([], f_s, dirichlet)

    
    def test_2d_td_solver(self):
        for n in [9]:
            dx, dy, dt = 1/(n+1), 1/(n+1), 1/800
            a = ddt.ddt(dt)
            b = td_d2dx.td_d2dx(dx, coefficient = lambda x, y: -1)
            c = td_d2dy.td_d2dy(dy, coefficient = lambda x, y: -1)
            expression = expr.diff_operator_expression([a, b, c])
            self.solver.diff_op_expression = expression
            time, result = self.solver.solve(n, n, 799)
            #print(time, result)
            assert len(time) == len(result)

if __name__ == '__main__':
    unittest.main()