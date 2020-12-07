import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../src/util/diff_operator_expression")
import diff_op_expression as expr
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../src/util/diff_operators/")
import impl.ddx as ddx, impl.ddy as ddy
import core.diff_op as diff_op
import unittest

class Test(unittest.TestCase):
    def setUp(self):
        self.ddx = ddx.ddx(0.1)
        self.expression = expr.diff_operator_expression()

    def test_append(self):
        self.expression.append(self.ddx)
    
    def test_multiple_append(self):
        self.expression.append([self.ddx, self.ddx])

if __name__ == '__main__':
    unittest.main()