from math import radians, sin, cos, tan

class Transformer:
    """
    * Transformation object for Ginzburg & Salmanova (1964) Globe Gore Projection
    * See Bugayevsky & Snyder 1995 p.222
    * NB: I have corrected this - the above reference version appears to be mis-specified
    """ 

    def __init__(self, R=6371000, k=2):
        """
        * constructor
        """
        self.R = R
        self.k = k

    def transform(self, lon, lat, direction='FORWARD'):
        """
        * Forward Transformation for Ginzburg & Salmanova (1964) Globe Gore Projection
        * See Bugayevsky & Snyder 1995 p.222
        * NB: I have corrected this - the above reference version appears to be mis-specified
        """ 
        # validate direction
        if (direction != "FORWARD"):
            raise NotImplementedError("Only the 'FORWARD' direction is currently implemented")
        
        # validate k
        if not 1 <= self.k <= 3:
            raise ValueError("parameter k should be in the range [1,3]")

        # convert to radians
        phi = radians(lat)
        lam = radians(lon)
        
        # deal with equator (avoid /0 error)
        if phi == 0: 
            x = lam
            y = 0
        
        # forward transform
        else:
            delta = lam * sin(phi) / self.k
            rho = self.k * (1 / tan(phi))
            x = rho * sin(delta)
            y = phi + rho * (1 - cos(delta))
        
        # scale & return
        return self.R * x, self.R * y