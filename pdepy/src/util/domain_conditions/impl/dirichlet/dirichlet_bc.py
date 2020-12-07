import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../..')
from domain_conditions.core.boundary_condition import bc # domain_conditions' core

class dirichlet_bc(bc.boundary_condition):
    def __init__(self, inDomain, onBoundary, getBoundaryValue, domain, getNearestPoint = (None, None)):
        super().__init__(inDomain, onBoundary, getBoundaryValue, getNearestPoint)
        self.domain = domain # this is the class domain