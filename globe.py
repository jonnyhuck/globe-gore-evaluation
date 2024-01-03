'''
* Evaluate projections for Globe gores
*
* @author jonnyhuck
'''

from pyproj import Geod
from numpy import arange
from pandas import DataFrame
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from pyproj import Geod, CRS, Transformer
from geopandas import read_file, GeoSeries
from distortion import evaluate_distortion, evaluate_gore_width, calibrate_metric, scale_factor

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

# construct a single gore of the specifided width at the Greenwich meridian
best_per_gore = []
for n_gores in [4, 6, 12, 24, 48]:
# for n_gores in [12]:
    gore_width = 360 / n_gores
    half_width = gore_width / 2
    gore = GeoSeries(Polygon([(-half_width, 90), (-half_width, -90), (half_width, -90), (half_width, 90)]), crs="EPSG:4326")

    # clip the land dataset ot the gore
    land = land.clip(gore)

    # get gore bounds
    minx, miny, maxx, maxy = gore.total_bounds

    # set geographical model against which distortion will be evaluated
    geo_str = "+proj=longlat +datum=WGS84 +no_defs"
    g = Geod(ellps='WGS84')

    # loop through each candidate projection
    results = []
    distributions = dict()
    width_plot = dict()
    for name, proj_str in crss.items():
        # print(name)
        
        # initialise a PyProj Transformer to transform coordinates
        transformer = Transformer.from_crs(CRS.from_proj4(geo_str), CRS.from_proj4(proj_str), always_xy=True)

        # evaluate gore width
        latitudes, w = evaluate_gore_width(transformer, gore_width, 0, 85, 5, 'mm', 400)

        # load into width plot
        width_plot[name] = {
            "x": w,
            "y": latitudes
        }

        # calculate distortion for current projection
        Ep, Es, Ea, Ka, p, s, a = evaluate_distortion(g, transformer, minx, miny+10, maxx, maxy-10, minr=10000, 
                                                      maxr=1000000, sample_number=10000, return_distributions=True)

        # load into distributions
        distributions[name] = {
            "Area": a,
            "Shape": s,
            "Distance": p,
        }

        # load into results
        results.append({
            "Name": name,
            "Area": Ea,
            "Shape": Es,
            "Distance": Ep,
        })
    

    # get max and min values
    max_vals = { "Area": 0, "Shape": 0, "Distance": 0 }
    min_vals = { "Area": 0, "Shape": 0, "Distance": 0 }
    for v in distributions.values():
        if max(v['Area']) > max_vals['Area']:
            max_vals['Area'] = max(v['Area'])
        if max(v['Shape']) > max_vals['Shape']:
            max_vals['Shape'] = max(v['Shape'])
        if max(v['Distance']) > max_vals['Distance']:
            max_vals['Distance'] = max(v['Distance'])

    # load into dataframe and run calibration on the results
    results = DataFrame(results)
    for id, row in results.iterrows():
        row['Area'] = scale_factor(calibrate_metric(max_vals['Area'], min_vals['Area'], Ea))
        row['Shape'] = calibrate_metric(max_vals['Shape'], min_vals['Shape'], Es)
        row['Distance'] = calibrate_metric(max_vals['Distance'], min_vals['Distance'], Ep)

    # calculate 'overall' distortion metric
    results['Total'] = (results['Distance'] + results['Shape'] + results['Area']) / 3
    
    # reporting...
    # print(f"Best for Distance: {results[results.Distance == results.Distance.min()]['Name'].iloc[0]}")
    # print(f"\t{results.sort_values('Distance')['Name'].to_list()}\n")
    # print(f"Best for Shape: {results[results.Shape == results.Shape.min()]['Name'].iloc[0]}")
    # print(f"\t{results.sort_values('Shape')['Name'].to_list()}\n")
    # print(f"Best for Area: {results[results.Area == results.Area.min()]['Name'].iloc[0]}")
    # print(f"\t{results.sort_values('Area')['Name'].to_list()}\n")
    # print(f"Best for Fit: {results[results.Width == results.Width.min()]['Name'].iloc[0]}")
    # print(f"\t{results.sort_values('Width')['Name'].to_list()}\n")
    print(f"Best overall for {n_gores} Gores: {results[results.Total == results.Total.min()]['Name'].iloc[0]}")
    print(f"\t{results.sort_values('Total')[['Name', 'Total', 'Area', 'Shape', 'Distance']]}\n")


    # plot distributions for distortion
    _, axs = plt.subplots(figsize=(17, 8), nrows=1, ncols=3)
    for ax, type in zip(axs, ["Area", "Shape", "Distance"]):
        # plt.ylabel(f"{type} Distortion")
        ax.set_title(f"{type} Distortion ({n_gores} Gores)")
        ax.boxplot(
            x=[distributions[n][type] for n in distributions.keys()], 
            flierprops=dict(markersize=6, markeredgecolor='#aaa', alpha=0.2),
            labels=list(crss.keys())
            )
        ax.yaxis.grid(True)
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
    plt.savefig(f"./out/{n_gores}_Gores.png", bbox_inches='tight')


    # plot lines for gore width error
    _, ax = plt.subplots(figsize=(5, 8), nrows=1, ncols=1)
    plt.axvline(x=0, color='black', lw=0.8)
    # ax.set_title(f"Gore 'Fit' error for a globe of 300 mm diameter")
    for name, data in width_plot.items():
        ax.plot(data['x'], data['y'], label=name)
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    plt.xlabel(f"Gore 'Fit' Error (mm)")
    plt.ylabel(f"Latitude")
    ax.yaxis.set_ticks(arange(0, 90, 10))
    plt.legend(loc="upper right")
    plt.savefig(f"./out/{n_gores}_gore_width.png", bbox_inches='tight')