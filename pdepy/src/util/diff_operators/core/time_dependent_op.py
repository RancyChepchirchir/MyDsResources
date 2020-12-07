from . import diff_op
class time_dependent_operator(diff_op.diff_operator):
    def __init__(self, implicit_stencil, explicit_stencil = None, coefficient = lambda x, y: 1):
        super().__init__(implicit_stencil, coefficient, True)
        self.explicit_stencil = explicit_stencil 
        # for 2d time-dependent equation. The first argument in tuple represents x coordinate offset,
        # while the second argument represents y coordinate offset

    def get_implicit_stencil(self):
        if self.stencil is None: raise NotImplementedError
        return self.stencil
    
    def get_explicit_stencil(self):
        if self.explicit_stencil is None: raise NotImplementedError
        return self.explicit_stencil