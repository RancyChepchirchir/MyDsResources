import numpy as np
from matplotlib.tri import Triangulation
from scipy.optimize import minimize

class poisson_solver:
    def __init__(self, f, domain_condition):
        self.f = f # if it's a time dependent problem, require f(x, t)
        self.domain_condition = domain_condition
        self.domain = domain_condition.domain
        self.tri = None
        self.points = None

    def solve(self, nx, ny):
        self.preprocess(nx, ny)

    def preprocess(self, nx, ny):
        def addPoint(x, y, isOnBoundary, X, Y):
            self.points.append((x, y, isOnBoundary))
            X.append(x)
            Y.append(y)
        dx, dy = self.domain.getDelta(nx, ny)
        X, Y = [], []
        self.points = []
        for j in range(0, ny+2): # grid index on y axis
            for i in range(0, nx+2): # grid index on x axis
                x, y = self._get_coord_by_offset(dx, dy, i, j) # the real (x,y) coordinate
                if self.domain_condition.inDomain(x, y):
                    addPoint(x, y, False, X, Y)
                elif self.domain_condition.onBoundary(x, y):
                    addPoint(x, y, True, X, Y)
        self.tri = Triangulation(X, Y).get_masked_triangles()


    def _get_coord_by_offset(self, dx, dy, x_grid_index, y_grid_index):
        x_offset, y_offset = x_grid_index, y_grid_index
        return self.domain.lower_left_coord + np.array([x_offset*dx, y_offset*dy])
