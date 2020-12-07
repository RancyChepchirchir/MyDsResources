class boundary_condition:
    def __init__(self, inDomain, onBoundary, getBoundaryValue, getNearestPoint = (None, None)):
        # inDomain is a function which takes a point, returns whether the point is in domain.
        self.inDomain = inDomain
        # onBoundary is a function which takes a point, returns whether the point is on the boundary of the domain.
        self.onBoundary = onBoundary
        # getBoundaryValue is a function which takes a point, and returns the boundary value on that point.
        self.getBoundaryValue = getBoundaryValue
        # get_nearest_point is a 2 element function tuple, 
        # where the first one takes (x, y) and gives the closest point on x axis,
        # and the second one gives the closest point on y axis.
        # used by irregular domain
        self.getNearestPoint = getNearestPoint