# Indicates that the domain can only be rectangle.
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../..')
from domain_conditions.core.boundary_condition import bc # domain_conditions' core
# from core.domain import domain as dm

class dirichlet_rectangular_bc(bc.boundary_condition):
    def __init__(self, inDomain, onBoundary, getBoundaryValue, domain, total_time, getNearestPoint = (None, None)):
        super().__init__(inDomain, onBoundary, getBoundaryValue, getNearestPoint)
        self.domain = domain # this is the class domain
        self.total_time = total_time