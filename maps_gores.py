'''
* Create comparison of 90 degree gores with approximated Tissot indicatrix
* 
* @author jonnyhuck
'''

from pyproj import Geod
import matplotlib.pyplot as plt
from geopandas import read_file, GeoSeries
from globe_gore_projection import Transformer as gTransformer
from shapely.geometry import Polygon, Point, LinearRing, LineString

# init globe gore transformer
gt = gTransformer()

def plot_globe_gore_projection(geometry):
    """
    * Recursive function to fool Geopandas into using a non-Proj transformation
    """
    if geometry.geom_type == 'Point':
        x, y = gt.transform(geometry.x, geometry.y)
        return Point(x, y)
    
    elif geometry.geom_type == 'LineString':
        # Transform each point in the LineString
        transformed_coords = [gt.transform(x, y) for x, y in geometry.coords]
        return LineString(transformed_coords)
    
    elif geometry.geom_type == 'LinearRing':
        # Transform each point in the LinearRing
        transformed_coords = [gt.transform(x, y) for x, y in geometry.coords]
        return LinearRing(transformed_coords)
    
    elif geometry.geom_type == 'Polygon':
        # Transform each ring in the Polygon (exterior and interior)
        exterior = [gt.transform(x, y) for x, y in geometry.exterior.coords]
        interiors = [[gt.transform(x, y) for x, y in interior.coords] for interior in geometry.interiors]
        return Polygon(exterior, interiors)
    
    else:
        raise TypeError(f"Geometry type '{geometry.geom_type}' is not supported")


def get_rhumb(startlong, startlat, endlong, endlat, nPoints):
    """
    * Calculates a rhumb line (required to draw boundary around gores)
    """

    # calculate distance between points
    g = Geod(f"+proj=latlong +datum=WGS84 +a=637100 +b=637100")

    # calculate line string along path with segments <= 1 km
    lonlats = g.npts(startlong, startlat, endlong, endlat, nPoints)

    # npts doesn't include start/end points, so prepend/append them
    lonlats.insert(0, (startlong, startlat))
    lonlats.append((endlong, endlat))
    
    # return the rhumb line
    return lonlats


# list of projections to evaluate
crss = {
    "Cassini": "+proj=cass +R=6371000",                    # equidistant
    "Polyconic": "+proj=poly +R=6371000",                  # compromise
    "Rectangular Polyconic": "+proj=rpoly +R=6371000",     # compromise
    "Ginzburg & Salmanova": "GLOBE",                       # equidistant
    "Sinusoidal": "+proj=sinu +R=6371000",                 # equal area
    "Transverse Mercator": "+proj=tmerc +R=6371000",       # conformal 
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
# TODO: convert
fig, axs = plt.subplots(figsize=(10, 10), nrows=2, ncols=3)
i = 0
crs_list = list(crss.items())
for ax_row in range(axs.shape[0]):
    for ax_col in range(axs.shape[1]):
        ax = axs[ax_row][ax_col]
        name, proj_str = crs_list[i]
        print(name)
        ax.axis('off')
        ax.set_title(f"{name}", pad= 20 if i==4 else 0 if i==2 else 10)
        
        # project using proj
        if proj_str != 'GLOBE':
            land.to_crs(proj_str).plot(ax=ax, facecolor='#aaa')
            graticule.to_crs(proj_str).plot(ax=ax, color='#ddd', lw=0.8)
            buffers.to_crs(proj_str).plot(ax=ax, facecolor='red', alpha=0.4)
            outline.to_crs(proj_str).plot(ax=ax, color='black')
        
        # fudge a custom projection
        else:
            land.geometry.apply(plot_globe_gore_projection).plot(ax=ax, facecolor='#aaa')
            graticule.geometry.apply(plot_globe_gore_projection).plot(ax=ax, color='#ddd', lw=0.8)
            buffers.geometry.apply(plot_globe_gore_projection).plot(ax=ax, facecolor='red', alpha=0.4)
            outline.geometry.apply(plot_globe_gore_projection).plot(ax=ax, color='black')
        i+=1

plt.subplots_adjust(wspace=0, hspace=0.2)
plt.savefig(f"./out/maps.png", bbox_inches='tight')