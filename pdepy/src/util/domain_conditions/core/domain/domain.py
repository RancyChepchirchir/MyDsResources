class domain:
    """
    Describes a crude domain where PDE exists. 
    Don't need to be too accurate as long as the real domain is a subset of this domain.
    This domain is a rectangle.
    """
    def __init__(self, lower_left_coord, upper_right_coord):
        # the lower left corner's coordinate of the rectangular domain. a two element tuple
        self.lower_left_coord = lower_left_coord
        # the upper right corner's coordinate of the rectangular domain. a two element tuple
        self.upper_right_coord = upper_right_coord

    def getDelta(self, nx, ny):
        x1, y1 = self.lower_left_coord
        x2, y2 = self.upper_right_coord
        return (abs(x2-x1)/(nx+1), abs(y2-y1)/(ny+1)) if ny != 1 else (abs(x2-x1)/(nx+1), None)