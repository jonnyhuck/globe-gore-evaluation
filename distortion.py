"""
* Various functions for evaluating map projection distortion
*
* If you run this file directly, it will print it will print out ideal proportional gore 
*    widths by latitude
*
* @author jonnyhuck
"""

from numpy import arange
from numpy.random import uniform
from shapely.geometry import Polygon
from math import hypot, sin, cos, radians, atan2, sqrt, pi


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


def spherical_haversine_distance(x1, y1, x2, y2, R = 6371000):
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


def evaluate_distortion(g, transformer, minx, miny, maxx, maxy, minr, maxr, sample_number, return_distributions=False):
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
    forward_azimuths = arange(0, 360, 22.5)
    n = len(forward_azimuths)

    # loop through the points
    planar_areas = []
    area_indices = []
    shape_indices = []
    distance_indices = []
    ellipsoidal_areas = []
    for x, y, r in zip(xs, ys, rs):

        # construct a circle around the centre point on the ellipsoid
        lons, lats = g.fwd([x]*n, [y]*n, forward_azimuths, [r]*n)[:2]

        # project the result, calculate area, append to the list
        e_coords = [ transformer.transform(lon, lat, direction='FORWARD') for lon, lat in zip(lons, lats) ]
        ellipsoidal_areas.append(Polygon(e_coords).area)

        # transform the centre point to the projected CRS
        px, py = transformer.transform(x, y, direction='FORWARD')

        # construct a circle around the projected point on a plane, calculate area, append to list
        p_coords = [ offset(px, py, r, az) for az in forward_azimuths ]
        planar_areas.append(Polygon(p_coords).area)

        # get radial distances frpm the centre to each of the 16 points on the circle
        ellipsoidal_radial_distances = [ hypot(px - ex, py - ey) for ex, ey in e_coords ]

        # get the sum of the distances, and the expected value for each distance
        total_radial_dist = sum(ellipsoidal_radial_distances)
        expected_distance = total_radial_dist / n

        # get the difference between the actual and expected radial distance for each 'spoke'
        shape_distortion = [ abs((expected_distance / total_radial_dist) - (d / total_radial_dist)) for d in ellipsoidal_radial_distances ]
        shape_indices.append(sum(shape_distortion))

    # calculate shape distortion
    Es = sum(shape_indices) / len(shape_indices)

    # calculate areal distortion (convert to Ka for each index, not all together at the end)
    diff_sum = 0
    for e, p in zip(ellipsoidal_areas, planar_areas):
        area_indices.append(abs(e - p) / abs(e + p))
        # diff_sum += scale_factor(abs(e - p) / abs(e + p))
        diff_sum += abs(e - p) / abs(e + p)
    Ea = 1 / sample_number * diff_sum


    ''' FINITE DISTANCE DISTORTION - as per Gosling & Symeonakis (2020) '''

    # loop once per sample required
    planar_distances = []
    ellipsoidal_distances = []
    for _ in range(sample_number):

        # get two random locations (x and y separately)
        xs = uniform(low=minx, high=maxx, size=2)
        ys = uniform(low=miny, high=maxy, size=2)

        # calculate the distance along the ellipsoid
        ellipsoidal_distances.append(g.line_length(xs, ys))

        # transform the coordinates
        origin = transformer.transform(xs[0], ys[0], direction='FORWARD')
        destination = transformer.transform(xs[1], ys[1], direction='FORWARD')

        # calculate the planar distance
        planar_distances.append(hypot(origin[0] - destination[0], origin[1] - destination[1]))

    # calculate distance distortion
    diff_sum = 0
    for e, p in zip(ellipsoidal_distances, planar_distances):
        distance_indices.append(abs(e - p) / abs (e + p))
        diff_sum += distance_indices[-1]
    Ep = 1 / sample_number * diff_sum

    # return all of the measures (and data if required)
    if return_distributions:
        return Ep, Es, Ea, distance_indices, shape_indices, area_indices
    else:
        return Ep, Es, Ea


def evaluate_gore_width(transformer, gore_width, min_lat=0, max_lat=90, interval=10, units="perc", globe_diameter=None):
    """
    * Compare the width of each gore at each interval of latitude to the same distance on the 
    *    sphere (not ellipsoid, as we are evaluating for the globe!)
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

        # get projected points and distance
        x1, y1 = transformer.transform(-half_width, lat, direction='FORWARD')
        x2, y2 = transformer.transform(half_width, lat, direction='FORWARD')
        projected_dist = hypot(x1-x2, y1-y2)

        # report the difference as a percentage of the gore
        if units == 'prop':
            # this is the difference as a proportion of the gore width
            distance_indices.append((projected_dist - spherical_dist) / spherical_dist)
        
        elif units == 'm':
            # this is the difference in m (in the world)
            distance_indices.append(projected_dist - spherical_dist)
        
        elif units == 'mm':
            # make sure we have a globe diameter
            if globe_diameter is None:
                print("units must be one of `prop`, `m`, `mm`")
                exit()

            # this is the difference in mm (on the printed globe)
            distance_indices.append((((projected_dist - spherical_dist) / spherical_dist) *  # difference as a proportion
                                    (2 * pi * globe_diameter / 2)) /                         # circumfrence of the globe
                                    (360 / gore_width)                                       # proportion of the globe copvered by the gore
                                    )

        else:
            print("units must be one of perc, m, mm")
            exit()

    # return total absolute error (and data to plot if required)
    return latitudes, distance_indices
    

def get_optimal_gore_width(interval=5):
    """
    * For use with flex projector to produce a 'perfect' gore projection
    * 
    * Multiply this number by the globe circumfrence to get the size on paper, or 
    *    by the gore width at the equator to get the actual width of the gores.
    """

    # get spherical distance across globe
    equator_dist = \
        spherical_haversine_distance(0, 0, 90, 0) + \
        spherical_haversine_distance(90, 0, 180, 0) + \
        spherical_haversine_distance(180, 0, -90, 0) + \
        spherical_haversine_distance(-90, 0, 0, 0)
    
    # for each latitude in interval
    for lat in list(range(0, 91, interval)):
        
        # get length of parallel
        spherical_dist = \
            spherical_haversine_distance(0, lat, 90, lat) + \
            spherical_haversine_distance(90, lat, 180, lat) + \
            spherical_haversine_distance(180, lat, -90, lat) + \
            spherical_haversine_distance(-90, lat, 0, lat)
    
        # report as proportion of equator dist
        print(f"{lat}: {spherical_dist / equator_dist}")

if __name__ == "__main__":
    """
    * If you run this file directly, it will print out ideal proportional gore 
    *   widths by latitude
    """
    get_optimal_gore_width()