import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/fdm/solver/")
import fdm_solver as fdm
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/diff_operator_expression")
import diff_op_expression as expr
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/")
import diff_operators.impl.ddx as ddx, diff_operators.impl.ddy as ddy, diff_operators.impl.laplacian2d as laplacian2d,\
    diff_operators.impl.ddt as ddt, diff_operators.impl.time_dependent_d2dx as td_d2dx
import diff_operators.core.diff_op as diff_op
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/core/domain")
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/impl/dirichlet")
import dirichlet_bc as dr
import domain as dm
import unittest
import numpy as np
import math

class Test(unittest.TestCase):
    def setUp(self):
        td_domain = dm.domain(np.array([0, 0]), np.array([1, 1]))
        td_inDomain = lambda x, y: 0 < x < 1 and 0 < y < 1
        td_onBoundary = lambda x, y: abs(x-1) < np.spacing(1) or abs(x) < np.spacing(1) \
            or abs(y) < np.spacing(1)
        def td_getBV(x, y):
            if abs(x) < np.spacing(1) or abs(x-1) < np.spacing(1):
                return 0
            else:
                #return math.sin(math.pi*x)
                #return 6*math.sin(math.pi*x) + 3*math.cos(2*math.pi*x) - 3
                return 6*x * math.sin(5*math.pi*x)
        f_s = lambda x, y: 0
        dirichlet = dr.dirichlet_bc(td_inDomain, td_onBoundary, td_getBV, td_domain)
        self.solver = fdm.fdm_solver([], f_s, dirichlet)

    def test_is_time_dependent(self):
        a = ddt.ddt(0.1)
        b = td_d2dx.td_d2dx(0.1)
        expression = expr.diff_operator_expression([a, b])
        assert expression.is_time_dependent() is True
        assert expression.all_ops_are_time_dependent() is True
        
    def test_all_ops_are_time_dependent(self):
        a = ddt.ddt(0.1)
        b = ddx.ddx(0.1)
        expression = expr.diff_operator_expression([a, b])
        assert expression.all_ops_are_time_dependent() is False
        assert expression.get_largest_coefficient() == 1
    
    def test_td_solver(self):
        for n in [50]:
            dx, dt = 1/(n+1), 1/(n*5+1)
            a = ddt.ddt(dt)
            b = td_d2dx.td_d2dx(dx, coefficient = lambda x, y: -1)
            expression = expr.diff_operator_expression([a, b])
            self.solver.diff_op_expression = expression
            #self.solver.solve(n, n)
            time, result = self.solver.solve(n, n*25)
            print(time, result)
            assert len(time) == len(result)

    def plot_time_dependent_function(self):
        # Uncomment the following code if GUI is supported
        # for i in range(len(a)):
        #     x1,x2,y1,y2 = plt.axis()
        #     plt.axis((x1,x2,0,6))
        #     plt.plot(np.linspace(0, 1, 41), a[i])
        #     plt.draw()
        #     plt.pause(0.001)
        #     plt.clf()
        pass

if __name__ == '__main__':
    unittest.main()