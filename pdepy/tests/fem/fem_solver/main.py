import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/fem/poisson_solver/")
import poisson
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/core/domain")
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/impl/dirichlet")
import dirichlet_bc as dr
import domain as dm
import unittest
import numpy as np
import math


class Test(unittest.TestCase):
    def setUp(self):
        domain = dm.domain(np.array([0, 0]), np.array([1, 1]))
        inDomain = lambda x, y: 0 < x < 1 and 0 < y < 1
        onBoundary = lambda x, y: abs(x) < np.spacing(1) or abs(x-1) < np.spacing(1) \
            or abs(y-1) < np.spacing(1) or abs(y) < np.spacing(1)
        getBoundaryValue = lambda x, y: 0
        diri = dr.dirichlet_bc(inDomain, onBoundary, getBoundaryValue, domain)
        self.solver = poisson.poisson_solver(lambda x, y: 32*x*(1-x)+32*y*(1-y), diri)
    
    def test_preprocess(self):
        self.solver.preprocess(3, 3)
        assert(len(self.solver.tri) == 32)
        assert(len(self.solver.points) == 25)

if __name__ == '__main__':
    unittest.main()