import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../util/diff_operators')
import core.time_dependent_op as td
from core.diff_op import Direction
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from math import ceil
class fdm_solver:
    def __init__(self, diff_op_expression, f, domain_condition):
        self.diff_op_expression = diff_op_expression
        self.f = f # if it's a time dependent problem, require f(x, t)
        self.domain_condition = domain_condition
        self.domain = domain_condition.domain
        self.index_to_grid = None
        self.grid_to_index = {}
        self.vector_len = 0
    
    def solve(self, nx, ny = 1, nt = None):
        """
        nx is the grid points number in x axis; ny is the symmetry of nx.
        The following grid has nx = 4, ny = 3
         _ _ _ _ _
        |_|_|_|_|_|
        |_|_|_|_|_|
        |_|_|_|_|_|
        |_|_|_|_|_|

        for 1d time-dependent problem, ny counts as nt, and nt will be none by default.
        for 2d time-dependent problem, ny and nt need to be set.
        """
        # For time dependent PDE, it requires that self.domain is exactly the domain of PDE
        if self.diff_op_expression.is_time_dependent():
            if not self.diff_op_expression.all_ops_are_time_dependent(): raise TypeError
            if nt is None: return self.time_dependent_implicit_solve(nx, ny)
            return self.time_dependent_explicit_solve(nx, ny, nt)
        else:
            return self.spatial_solve(nx, ny)
    
    def spatial_solve(self, nx, ny):
        # This function handles the case where all the operators are not time-dependent.
        self.preprocess(nx, ny)
        dx, dy = self.domain.getDelta(nx, ny)
        A = np.zeros([self.vector_len, self.vector_len])
        fv, u = np.zeros(self.vector_len), np.zeros(self.vector_len)
        for index in range(self.vector_len):
            x, y = self.index_to_grid[index]
            fv[index] = self.f(*self._get_coord_by_offset(dx, dy, x, y))
            x_coord, y_coord = self._get_coord_by_offset(dx, dy, x, y)
            for op in self.diff_op_expression:
                history_u, history_A = [], [] # record operations for each operator in order to revert back
                for node in op.stencil:
                    coord, coeff = node
                    x_offset, y_offset = coord
                    cur_x, cur_y = x + x_offset, y + y_offset
                    cur_x_coord, cur_y_coord = self._get_coord_by_offset(dx, dy, cur_x, cur_y)
                    if self.domain_condition.onBoundary(cur_x_coord, cur_y_coord):
                        bv = self.domain_condition.getBoundaryValue(cur_x_coord, cur_y_coord)
                        u[index] -= op.coefficient(cur_x_coord, cur_y_coord) * coeff * bv
                        history_u.append((index, op.coefficient(cur_x_coord, cur_y_coord) * coeff * bv))
                    elif not self.domain_condition.inDomain(cur_x_coord, cur_y_coord):
                        if not self._has_getNearestPoint(ny): raise NotImplementedError
                        self._op_revert_back(history_A, history_u, A, u)
                        if x_offset != 0:
                            direction = Direction.POSITIVE if x_offset > 0 else Direction.NEGATIVE
                            bx_coord, by_coord = self.domain_condition.getNearestPoint[0](cur_x_coord, cur_y_coord)
                            tau = abs(bx_coord - x_coord) / dx
                            irregular_stencil = op.getIrregularStencil(dx, tau, direction)
                            for irregular_node in irregular_stencil:
                                x_offset, y_offset = irregular_node.offset
                                cur_x, cur_y = x + x_offset, y + y_offset
                                cur_x_coord, cur_y_coord = self._get_coord_by_offset(dx, dy, cur_x, cur_y)
                                if irregular_node.flawed:
                                    bv = self.domain_condition.getBoundaryValue(cur_x_coord, cur_y_coord)
                                    u[index] -= op.coefficient(cur_x_coord, cur_y_coord) * irregular_node.coefficient * bv
                                else:
                                    A[index, self.grid_to_index[cur_x, cur_y]] += op.coefficient(cur_x_coord, cur_y_coord) * irregular_node.coefficient
                        elif y_offset != 0:
                            direction = Direction.POSITIVE if y_offset > 0 else Direction.NEGATIVE
                            bx_coord, by_coord = self.domain_condition.getNearestPoint[1](cur_x_coord, cur_y_coord)
                            tau = abs(by_coord - y_coord) / dy
                            irregular_stencil = op.getIrregularStencil(dy, tau, direction)
                            for irregular_node in irregular_stencil:
                                x_offset, y_offset = irregular_node.offset
                                cur_x, cur_y = x + x_offset, y + y_offset
                                cur_x_coord, cur_y_coord = self._get_coord_by_offset(dx, dy, cur_x, cur_y)
                                if irregular_node.flawed:
                                    bv = self.domain_condition.getBoundaryValue(cur_x_coord, cur_y_coord)
                                    u[index] -= op.coefficient(cur_x_coord, cur_y_coord) * irregular_node.coefficient * bv
                                else:
                                    A[index, self.grid_to_index[cur_x, cur_y]] += op.coefficient(cur_x_coord, cur_y_coord) * irregular_node.coefficient
                        else:
                            raise TypeError
                        break
                    else:
                        A[index, self.grid_to_index[cur_x, cur_y]] += op.coefficient(cur_x_coord, cur_y_coord) * coeff
                        history_A.append((index, self.grid_to_index[cur_x, cur_y], -op.coefficient(cur_x_coord, cur_y_coord) * coeff))
        return spsolve(csr_matrix(A), fv + u)
        
    def time_dependent_implicit_solve(self, nx, nt):
        dx, dt = self.domain.getDelta(nx, nt)
        result = self._td_get_initial_value(nx)
        time = [0]
        for j in range(1, nt+2):
            A = np.zeros([nx, nx])
            fv, u = np.zeros(nx), np.zeros(nx)
            for i in range(1, nx+1): # from 1 to nx
                fv[i-1] = self.f(*self._get_coord_by_offset(dx, dt, i-1, j-1))
                for op in self.diff_op_expression:
                    for node in op.stencil:
                        coord, coeff = node
                        x_offset, t_offset = coord
                        cur_x, cur_t = i + x_offset, j + t_offset
                        cur_x_coord, cur_t_coord = self._get_coord_by_offset(dx, dt, cur_x-1, cur_t-1)
                        if self.domain_condition.onBoundary(cur_x_coord, cur_t_coord):
                            bv = self.domain_condition.getBoundaryValue(cur_x_coord, cur_t_coord)
                            u[i-1] -= op.coefficient(cur_x_coord, cur_t_coord) * coeff * bv
                        elif t_offset < 0: # the value has already been computed in the previous computations
                            u[i-1] -= op.coefficient(cur_x_coord, cur_t_coord) * coeff * result[cur_t, cur_x]
                        else:
                            A[i-1, cur_x-1] += op.coefficient(cur_x_coord, cur_t_coord) * coeff
            row = spsolve(csr_matrix(A), fv + u)
            x1, y = self._get_coord_by_offset(dx, dt, -1, j-1)
            x2 = self._get_coord_by_offset(dx, dt, nx, j-1)[0]
            f1, f2 = self.domain_condition.getBoundaryValue(x1, y), self.domain_condition.getBoundaryValue(x2, y)
            row = np.concatenate([[f1], row, [f2]])
            result = np.vstack([result, row])
            time.append(j*dt)
        return time, result

    def time_dependent_explicit_solve(self, nx, ny, nt):
        if self.domain_condition.total_time is None: raise NotImplementedError
        dx, dy = self.domain.getDelta(nx, ny)
        max_coefficient = self.diff_op_expression.get_largest_coefficient()
        dt = self.domain_condition.total_time / (nt+1)
        if dt > (dx*dy)**2/(2*max_coefficient*(dx**2+dy**2)):
            recommended_nt = ceil(self.domain_condition.total_time/((dx*dy)**2/(2*max_coefficient*(dx**2+dy**2))))
            sys.stderr.write(f"Since the program uses explicit method to solve 2D time-dependent PDE, it's recommended to set nt to be larger than {recommended_nt}\n")
        last = self._td_get_2d_initial_value(nx, ny)
        result, time = [last], [0]
        lower_left_x, lower_left_y = self.domain.lower_left_coord
        upper_right_x, upper_right_y = self.domain.upper_right_coord
        X, Y = np.linspace(lower_left_x, upper_right_x, nx+2), np.linspace(lower_left_y, upper_right_y, ny+2)
        for j in range(1, nt+2):
            current = np.zeros([ny+2, nx+2])
            for y_grid_index in range(0, ny+2):
                for x_grid_index in range(0, nx+2):
                    x, y = X[x_grid_index], Y[y_grid_index]
                    if self.domain_condition.onBoundary(x, y):
                        current[y_grid_index][x_grid_index] = self.domain_condition.getBoundaryValue(x, y)
                    else:
                        left_coeff, right_const = 0, 0
                        for op in self.diff_op_expression:
                            for node in op.get_explicit_stencil():
                                coord, coeff = node
                                cur_x_index, cur_y_index = x_grid_index + coord[0], y_grid_index + coord[1]
                                op_type = coord[2]
                                if op_type == -1:
                                    right_const -= op.coefficient(X[cur_x_index], Y[cur_y_index]) * coeff * last[cur_y_index][cur_x_index]
                                elif op_type == 0:
                                    left_coeff += op.coefficient(X[cur_x_index], Y[cur_y_index]) * coeff
                                else:
                                    raise NotImplementedError
                        right_const += self.f(x, y)
                        current[y_grid_index][x_grid_index] = right_const / left_coeff
            result.append(current)
            last = current
            time.append(j*dt)
        return time, result

    def _op_revert_back(self, history_A, history_u, A, u):
        for c_x, c_y, v in history_A:
            A[c_x, c_y] += v
        for index, v in history_u:
            u[index] += v

    def _td_get_2d_initial_value(self, nx, ny):
        lower_left_x, lower_left_y = self.domain.lower_left_coord
        upper_right_x, upper_right_y = self.domain.upper_right_coord
        x, y = np.linspace(lower_left_x, upper_right_x, nx+2), np.linspace(lower_left_y, upper_right_y, ny+2)
        result = np.zeros([ny+2, nx+2])
        for y_grid_index in range(0, ny+2):
            for x_grid_index in range(0, nx+2):
                result[y_grid_index][x_grid_index] = self.domain_condition.getBoundaryValue(x[x_grid_index], y[y_grid_index])
        return result

    def _td_get_initial_value(self, nx):
        """
        Compute the initial value of the time dependent pde problem
        """
        lower_left_x, lower_left_y = self.domain.lower_left_coord
        upper_right_x = self.domain.upper_right_coord[0]
        x, y = np.linspace(lower_left_x, upper_right_x, nx+2), np.empty(nx+2)
        y.fill(lower_left_y)
        result = np.zeros(nx+2)
        for i in range(nx+2):
            result[i] = self.domain_condition.getBoundaryValue(x[i], y[i])
        return result.reshape(1, nx+2)

    def preprocess(self, nx, ny = 1):
        dx, dy = self.domain.getDelta(nx, ny)
        self.index_to_grid = np.zeros(nx*ny, dtype=(int, 2))
        self.vector_len = 0
        for j in range(0, ny): # grid index on y axis
            for i in range(0, nx): # grid index on x axis
                x, y = self._get_coord_by_offset(dx, dy, i, j) # the real (x,y) coordinate
                if self.domain_condition.inDomain(x, y):
                    self.index_to_grid[self.vector_len] = (i, j)
                    self.grid_to_index[(i, j)] = self.vector_len
                    self.vector_len += 1

    def _get_coord_by_offset(self, dx, dy, x_grid_index, y_grid_index):
        x_offset, y_offset = x_grid_index + 1, y_grid_index + 1
        return self.domain.lower_left_coord + np.array([x_offset*dx, y_offset*dy])

    def _has_getNearestPoint(self, ny):
        if self.domain_condition.getNearestPoint[0] is None:
            return False
        elif self.domain_condition.getNearestPoint[1] is None and ny != 1:
            return False
        return True