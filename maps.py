'''
* Create comparison of 90 degree gores with approximated Tissot indicatrix
* 
* @author jonnyhuck
'''

from pyproj import Geod
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, LinearRing
from geopandas import read_file, GeoSeries


def get_rhumb(startlong, startlat, endlong, endlat, nPoints):
    """
    * Calculates a rhumb line (required to draw boundary around gores)
    """

    # calculate distance between points
    g = Geod(ellps='WGS84')

    # calculate line string along path with segments <= 1 km
    lonlats = g.npts(startlong, startlat, endlong, endlat, nPoints)

    # npts doesn't include start/end points, so prepend/append them
    lonlats.insert(0, (startlong, startlat))
    lonlats.append((endlong, endlat))
    
    # return the rhumb line
    return lonlats


# list of projections to evaluate
crss = {
    "Sinusoidal": "+proj=sinu +R=6371000",                          # equal area
    "Cassini": "+proj=cass +R=6371000",                             # equidistant - transverse of Plate Carree
    "Transverse Mercator": "+proj=tmerc +R=6371000",                # conformal 
    "Rectangular Polyconic": "+proj=rpoly +R=6371000",              # compromise
    "American Polyconic": "+proj=poly +R=6371000",                  # compromise
    }

# open landmass data
land = read_file("./data/ne_110m_land.shp")
graticule = read_file("./data/ne_110m_graticules_30.shp")

# construct a single gore of the specifided width at the Greenwich meridian
n_gores = 4
gore_width = 360 / n_gores
half_width = gore_width / 2
minx, miny, maxx, maxy = [-half_width, -90, half_width, 90.0]
lr = LinearRing(
    get_rhumb(maxx, maxy, maxx, 0, 500) +    # north pole to equator down east side
    get_rhumb(maxx, 0, maxx, miny, 500) +    # equator to south pole down east side
    get_rhumb(minx, miny, minx, 0, 500) +    # south pole to equator up west side
    get_rhumb(minx, 0, minx, maxy, 500)      # equator to north pole up west side
    )
outline = GeoSeries(lr, crs="EPSG:4326")
gore = GeoSeries(Polygon(lr), crs="EPSG:4326")

# create approximated tissot indicatrix (geodesic circles)
buffers = GeoSeries([Point((lng, lat)).buffer(8) for lng in range (int(minx), 
    int(maxx)+10, 22) for lat in range (int(miny), int(maxy)+10, 22)], crs="EPSG:4326")

# clip the land dataset at the gore
land = land.clip(gore)
graticule = graticule.clip(gore)
buffers = buffers.clip(gore)

# plot figure
_, axs = plt.subplots(figsize=(20, 8), nrows=1, ncols=5)
i = 0
crs_list = list(crss.items())
for ax in axs:
    name, proj_str = crs_list[i]
    print(name)
    ax.axis('off')
    if name != 'Transverse Mercator':
        ax.set_title(f"{name}")
    else:
        ax.set_title(f"{name}\n\n")
    land.to_crs(proj_str).plot(ax=ax, facecolor='#aaa')
    graticule.to_crs(proj_str).plot(ax=ax, color='#ddd', lw=0.8)
    buffers.to_crs(proj_str).plot(ax=ax, facecolor='red', alpha=0.4)
    outline.to_crs(proj_str).plot(ax=ax, color='black')
    i+=1
plt.savefig(f"./out/maps.png", bbox_inches='tight')