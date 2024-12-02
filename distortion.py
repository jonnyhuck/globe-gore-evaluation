"""
* Various functions for evaluating map projection distortion for SPHERICAL models
*
* If you run this file directly, it will print it will print out ideal proportional gore 
*    widths by latitude
*
* @author jonnyhuck
"""
from numpy import arange
from numpy.random import uniform
from shapely.geometry import Polygon
from math import hypot, sin, cos, radians, atan2, sqrt, pi, degrees, asin

def spherical_forward(lon1, lat1, bearing, distance, R=6371000):
    """
    * Calculate the destination point given an initial point, a bearing, and a distance on a sphere
    """
    # Convert latitude, longitude, and bearing to radians
    lon1 = radians(lon1)
    lat1 = radians(lat1)
    bearing = radians(bearing)

    # Angular distance in radians
    angular_distance = distance / R

    # Calculate the destination latitude
    lat2 = asin(sin(lat1) * cos(angular_distance) + cos(lat1) * sin(angular_distance) * cos(bearing))

    # Calculate the destination longitude
    lon2 = lon1 + atan2(sin(bearing) * sin(angular_distance) * cos(lat1), cos(angular_distance) - sin(lat1) * sin(lat2))

    # Convert to degrees, ensure correct range and return
    return (degrees(lon2)+ 180) % 360 - 180, degrees(lat2)


def sfwd(lon1s, lat1s, bearings, distances, R=6371000):
    """ 
    * Convenience version of the above to match functionality of g.fwd() 
    """
    out = [spherical_forward(lon1, lat1, bearing, distance, R) for lon1, lat1, bearing, distance in zip(lon1s, lat1s, bearings, distances)]
    lons = [o[0] for o in out]
    lats = [o[1] for o in out]
    return lons, lats


def offset(x, y, distance, direction):
    """
    * Offset a location by a given distance and direction on a plane
    """
    x2 = x + cos(radians(direction)) * distance
    y2 = y + sin(radians(direction)) * distance
    return (x2, y2)


def scale_factor(Ea):
    """
    * Convert Ea to scale factor, as per Canters et al. (2005)
    """
    return (1 + Ea) / (1 - Ea)


def calibrate_metric(E_max, E_min, E):
    """
    * Calibrate the distortion metrics to make them comparable, as per Canters et al. (2005)
    """
    return 1 / (E_max - E_min) * (E - E_min)


def spherical_haversine_distance(x1, y1, x2, y2, R=6371000):
    """
    Calculate the 'as the crow flies' distance between two locations along a sphere using
        the Haversine equation.
    """
    # convert to radians (as this is what the math functions expect)
    a1 = radians(y1)
    b1 = radians(x1)
    a2 = radians(y2)
    b2 = radians(x2)

    # half sine
    sa = sin((a2 - a1)*0.5)
    sb = sin((b2 - b1)*0.5)

    # the square of half the chord length between the points
    a = sa * sa + cos(a1) * cos(a2) * sb * sb

    # angular distance in degrees
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return c * float(R)


def sdist(lons, lats, R=6371000):
    """
    * Convenience version of the above to match functionality of g.line_length() 
    """
    return spherical_haversine_distance(lons[0], lats[0], lons[1], lats[1], R)


def evaluate_distortion(transformer, minx, miny, maxx, maxy, minr, maxr, sample_number, vertices=16, return_distributions=False):
    """
    * Calculate a selection of distortion measures, based on Canters et al. (2005)
    *  and Gosling & Symeonakis (2020)
    """

    ''' FINITE AREAL AND SHAPE DISTORTION - based on Canters et al. (2005) '''

    # calculate the required number of random locations (x and y separately) plus radius
    xs = uniform(low=minx, high=maxx, size=sample_number)
    ys = uniform(low=miny, high=maxy, size=sample_number)
    rs = uniform(low=minr, high=maxr, size=sample_number)

    # offset distances
    forward_azimuths = arange(0, 360, 360/vertices)

    # loop through the points
    area_indices = []
    shape_indices = []
    distance_indices = []
    for x, y, r in zip(xs, ys, rs):

        # construct a circle around the centre point on the sphere
        lons, lats = sfwd([x]*vertices, [y]*vertices, forward_azimuths, [r]*vertices)[:2]

        # project the result, calculate area, append to the list
        s_coords = [ transformer.transform(lon, lat, direction='FORWARD') for lon, lat in zip(lons, lats) ]
        spherical_area = Polygon(s_coords).area

        # transform the centre point to the projected CRS
        px, py = transformer.transform(x, y, direction='FORWARD')

        # construct a circle around the projected point on a plane, calculate area, append to list
        planar_area = Polygon([ offset(px, py, r, az) for az in forward_azimuths ]).area

        # get area index
        area_indices.append(abs(spherical_area - planar_area) / abs(spherical_area + planar_area))

        # get radial distances from the centre to each of the 16 points on the circle
        spherical_radial_distances = [ hypot(px - sx, py - sy) for sx, sy in s_coords ]

        # get the absolute proportional difference between the expected and actual radial distance for each 'spoke'
        shape_distortion = [abs((1 / vertices) - (d / sum(spherical_radial_distances))) for d in spherical_radial_distances]
        shape_indices.append(sum(shape_distortion))

    # calculate areal & shape distortion
    Ea = 1 / sample_number * sum(area_indices)      # as per equation
    Es = sum(shape_indices) / len(shape_indices)    # mean value

    
    ''' FINITE DISTANCE DISTORTION - as per Gosling & Symeonakis (2020) '''

    # loop once per sample required
    for _ in range(sample_number):

        # get two random locations (x and y separately)
        xs = uniform(low=minx, high=maxx, size=2)
        ys = uniform(low=miny, high=maxy, size=2)

        # calculate the distance along the ellipsoid
        spherical_distance = sdist(xs, ys)

        # transform the coordinates
        origin = transformer.transform(xs[0], ys[0], direction='FORWARD')
        destination = transformer.transform(xs[1], ys[1], direction='FORWARD')

        # calculate the planar distance
        planar_distance = hypot(origin[0] - destination[0], origin[1] - destination[1])
        
        # calculate distance index 
        distance_indices.append(abs(spherical_distance - planar_distance) / abs (spherical_distance + planar_distance))

    # calculate distance distortion
    Ep = 1 / sample_number * sum(distance_indices)

    # get min max for calibration if needed and return
    if return_distributions:
        return Ep, Es, Ea, distance_indices, shape_indices, area_indices
    
    # otherwise just return indices
    else:
        return Ep, Es, Ea


def evaluate_gore_fit(transformer, gore_width, sample_number, min_lat=0, max_lat=89, return_distributions=False):
    """
    * Provides an evaluation of gore fit to the sphere that is analogous to the measures used in 
    *   `evaluate_distortion()`. This provides a single value, as opposed to `evaluate_gore_width()`, 
    *   which evaluates by latitude.
    """
     # get the distance from the central meridian to the gore edge
    half_width = gore_width / 2

    # loop through each line of latitude
    distance_indices = []
    for lat in uniform(low=min_lat, high=max_lat, size=sample_number):

        # get spherical distance across gore
        spherical_distance = spherical_haversine_distance(-half_width, lat, half_width, lat)

        # get projected distance
        x1, y1 = transformer.transform(-half_width, lat, direction='FORWARD')
        x2, y2 = transformer.transform(half_width, lat, direction='FORWARD')
        projected_distance = hypot(x1-x2, y1-y2)

        # get distance index
        distance_indices.append(abs(spherical_distance - projected_distance) / abs (spherical_distance + projected_distance))
    
    # get fit metric
    Ef = 1 / sample_number * sum(distance_indices)

    # return total absolute error (and data to plot if required)
    if return_distributions:
        return Ef, distance_indices
    else:
        return Ef


def evaluate_gore_width(transformer, gore_width, min_lat=0, max_lat=90, interval=10, units="prop", globe_diameter=None):
    """
    * Compare the width of each gore at each interval of latitude to the same distance across the sphere
    * Calculations are made for each 
    """

    # get the distance from the central meridian to the gore edge
    half_width = gore_width / 2

    # loop through each line of latitude
    latitudes = []
    distance_indices = []
    for lat in range(min_lat, max_lat+1, interval):

        # store latitude 
        latitudes.append(lat)

        # get spherical distance across gore
        spherical_dist = spherical_haversine_distance(-half_width, lat, half_width, lat)

        # get projected distance
        x1, y1 = transformer.transform(-half_width, lat, direction='FORWARD')
        x2, y2 = transformer.transform(half_width, lat, direction='FORWARD')
        projected_dist = hypot(x1-x2, y1-y2)

        # report the difference as a percentage of the gore
        if units == 'prop':
            # this is the difference as a proportion of the gore width
            distance_indices.append((projected_dist - spherical_dist) / spherical_dist)
        
        elif units == 'm':
            # this is the difference in m (equivelant distance 'in the world')
            distance_indices.append(projected_dist - spherical_dist)
        
        elif units == 'mm':
            # make sure we have a globe diameter
            if globe_diameter is None:
                print("please provide a globe diameter!")
                exit()

            # this is the difference in mm (on the printed globe)
            distance_indices.append((((projected_dist - spherical_dist) / spherical_dist) *  # difference as a proportion
                                    (2 * pi * globe_diameter / 2)) /                         # circumfrence of the globe
                                    (360 / gore_width)                                       # proportion of the globe covered by the gore
                                    )
        else:
            print("units must be one of prop, m, mm")
            exit()

    # return total absolute error (and data to plot if required)
    return latitudes, distance_indices
    

# evaluate globe projection, just for testing
if __name__ == "__main__":
    from globe_gore_projection import Transformer as gTransformer
    transformer = gTransformer(6371000, 2)
    lats, errors = evaluate_gore_width(transformer, 360/12, max_lat=85, interval=5, units='mm', globe_diameter=500)
    for l, w in zip(lats, errors):
        print(l, w)