'''
* Evaluate projections for Globe gores
*
* This version outputs the box and whisker plots as png files and the individual 
*   scoring tables to the terminal
*
* @author jonnyhuck
'''
from numpy import arange
from pandas import DataFrame
from geopandas import GeoSeries
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
from shapely.geometry import Polygon
from globe_gore_projection import Transformer as gTransformer
from distortion import evaluate_distortion, evaluate_gore_width, evaluate_gore_fit, calibrate_metric

# earth radius
RADIUS = 6371000

# list of projections to evaluate
crss = {
    "Cassini":                  f"+proj=cass +R={RADIUS} +datum=WGS84",   # equidistant
    "Polyconic":                f"+proj=poly +R={RADIUS} +datum=WGS84",   # compromise
    "Rectangular Polyconic":    f"+proj=rpoly +R={RADIUS} +datum=WGS84",  # compromise
    "Ginzburg & Salmanova":     f"GLOBE",                                 # equidistant 
    "Sinusoidal":               f"+proj=sinu +R={RADIUS} +datum=WGS84",   # equal area
    "Transverse Mercator":      f"+proj=tmerc +R={RADIUS} +datum=WGS84",  # conformal 
    }

# weights for weighted average
weights = {
    "area": 0.25,
    "shape": 0.25,
    "distance": 0,
    "fit": 0.5,
}

# results dictionary
overall_distortion = dict()

# loop through all diameters
for diameter in [220, 360, 500, 650, 800]:
# for diameter in [500]:
    print(f"{diameter}mm globe:")

    # construct results dictionary
    overall_distortion[diameter] = dict()
    for c in crss.keys():
        overall_distortion[diameter][c] = dict()
        overall_distortion[diameter][c]['Distortion'] = []
        overall_distortion[diameter][c]['Total'] = []

    # evaluate all desired gore numbers
    for n_gores in [6, 12, 24]:
    # for n_gores in [12]:

        # set geographical model against which distortion will be evaluated
        geo_str = f"+proj=latlong +datum=WGS84 +a={RADIUS} +b={RADIUS}"

        # construct a single gore of the specifided width at the Greenwich meridian
        gore_width = 360 / n_gores
        half_width = gore_width / 2
        gore = GeoSeries(Polygon([(-half_width, 90), (-half_width, -90), (half_width, -90), (half_width, 90)]), crs=geo_str)

        # get gore bounds
        minx, miny, maxx, maxy = gore.total_bounds

        # loop through each candidate projection
        results = []
        distributions = dict()
        width_plot = dict()
        for name, proj_str in crss.items():
            
            # initialise a PyProj Transformer to transform coordinates
            if proj_str[:5] == "GLOBE":
                transformer = gTransformer(RADIUS, 2)
            else:
                transformer = Transformer.from_crs(CRS.from_proj4(geo_str), CRS.from_proj4(proj_str), always_xy=True)

            # evaluate gore width
            latitudes, w = evaluate_gore_width(transformer, gore_width, 0, 85, 5, 'mm', diameter)
            Ef, f = evaluate_gore_fit(transformer, gore_width, 10000, 0, 85, True)

            # load into width plot
            width_plot[name] = {
                "x": w,
                "y": latitudes
            }

            # calculate distortion for current projection
            Ep, Es, Ea, p, s, a = evaluate_distortion(transformer, minx, miny+10, maxx, maxy-10, minr=100000, 
                                                        maxr=1000000, sample_number=10000, return_distributions=True)

            # load into distributions
            distributions[name] = {
                "Area": a,
                "Shape": s,
                "Distance": p,
                "Fit": f,
            }

            # load into results
            results.append({
                "Name": name,
                "Area": Ea,
                "Shape": Es,
                "Distance": Ep,
                "Fit": Ef
            })

        # get max index returned for each projection/measure 
        max_vals = { "Area": -float('inf'), "Shape": -float('inf'), "Distance": -float('inf'), "Fit": -float('inf') }
        for v in distributions.values():
            max_vals['Area'] = max(max_vals['Area'], max(v['Area']))
            max_vals['Shape'] = max(max_vals['Shape'], max(v['Shape']))
            max_vals['Distance'] = max(max_vals['Distance'], max(v['Distance']))
            max_vals['Fit'] = max(max_vals['Fit'], max(v['Fit']))

        # calibrate results between 0 and max observed value
        for r in results:
            r['Area'] = calibrate_metric(max_vals['Area'], 0, r['Area'])
            r['Shape'] = calibrate_metric(max_vals['Shape'], 0, r['Shape'])
            r['Distance'] = calibrate_metric(max_vals['Distance'], 0, r['Distance'])
            r['Fit'] = calibrate_metric(max_vals['Fit'], 0, r['Fit'])

        # convert to dataframe and calculate 'overall' metric as weighted mean
        results = DataFrame(results)
        results['Total'] = (
                results['Distance'] * weights['distance'] + 
                results['Shape'] * weights['shape'] + 
                results['Area'] * weights['area'] +
                results['Fit'] * weights['fit']
            )
        
        # load into dictionary
        for id, row in results.iterrows():
            overall_distortion[diameter][row['Name']]['Total'].append(row.loc['Total'])
        
        # reporting...
        print(f"{results}\n")
        print(f"{n_gores} gores fit: {results[results.Fit == results.Fit.min()]['Name'].iloc[0]}; overall score: {results[results.Total == results.Total.min()]['Name'].iloc[0]}\n\n")

        # plot distributions for distortion
        _, axs = plt.subplots(figsize=(22, 8), nrows=1, ncols=4)
        for ax, type in zip(axs, ["Area", "Shape", "Distance", "Fit"]):
            ax.set_title(f"{type} {'Distortion ' if type != 'Fit' else ''}({n_gores} Gores, {diameter} mm Sphere)")
            ax.boxplot(
                x=[distributions[n][type] for n in distributions.keys()], 
                flierprops=dict(markersize=6, markeredgecolor='#aaa', alpha=0.2),
                labels=list(crss.keys())
                )
            ax.yaxis.grid(True)
            ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
        plt.savefig(f"./out/{diameter}mm_{n_gores}_gores.png", bbox_inches='tight')

        # plot lines for gore width error
        _, ax = plt.subplots(figsize=(5, 8), nrows=1, ncols=1)
        plt.axvline(x=0, color='black', lw=0.8)
        for name, data in width_plot.items():
            ax.plot(data['x'], data['y'], label=name)
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)
        plt.xlabel(f"Gore 'Fit' Error (mm)")
        plt.ylabel(f"Latitude")
        ax.yaxis.set_ticks(arange(0, 90, 10))
        plt.legend(loc="upper right")
        plt.savefig(f"./out/{diameter}mm_{n_gores}_gore_width.png", bbox_inches='tight')
        plt.close('all')