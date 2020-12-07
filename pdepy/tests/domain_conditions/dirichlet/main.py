import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/impl/dirichlet")
import dirichlet_bc as dr
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/core/domain")
import domain as dm
import numpy as np
import unittest

class Test(unittest.TestCase):
    def setUp(self):
        self.inDomain = lambda x, y: abs(x) < 1 and abs(y) < 1
        self.onBoundary = lambda x, y: abs(x-1) < np.spacing(1) or abs(x+1) < np.spacing(1) \
            or abs(y-1) < np.spacing(1) or abs(y+1) < np.spacing(1)
        self.getBoundaryValue = lambda x, y: 1
        self.domain = dm.domain(np.array([-1, -1]), np.array([1, 1]))

    def test_init(self):
        a = dr.dirichlet_bc(self.inDomain, self.onBoundary, self.getBoundaryValue, self.domain)
        assert (a.domain.lower_left_coord == np.array([-1, -1])).all()

if __name__ == '__main__':
    unittest.main()