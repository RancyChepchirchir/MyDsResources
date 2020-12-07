import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../../src/util/domain_conditions/core/domain")
import domain as dm
import unittest
import numpy as np

class Test(unittest.TestCase):
    def setUp(self):
        self.domain = dm.domain(np.array([-1, -1]), np.array([1, 1]))

    def test_init(self):
        delta_x, delta_y = self.domain.getDelta(9, 19)
        assert abs(delta_x - 0.2) < np.spacing(1)
        assert abs(delta_y - 0.1) < np.spacing(1)
    
    def test_1d(self):
        delta_x, delta_y = self.domain.getDelta(9, 1)
        assert abs(delta_x - 0.2) < np.spacing(1)
        assert delta_y == None

if __name__ == '__main__':
    unittest.main()