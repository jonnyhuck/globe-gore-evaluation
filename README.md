# globe-gore-evaluation

Evaluating the distortion and fit of candidate projections for globe gores

Following my [earlier interest in globes](https://jonnyhuckblog.wordpress.com/2016/06/29/globemaking-for-beginners/) and (now outdated) [code to produce one](https://github.com/jonnyhuck/GlobeMaker) - this repository provides the code and results of the first rigorous analysis pf the properties of a range of projections for use in globemaking.

* `distortion.py` contains library code for evaluating the projections
* `globe.py` contains the evaluation itself
* `maps.py` simply produces a figure comparing the projections under consideration

This code allow you to evaluate the distortion properties (shown here for a 12-gore globe):

![distortion properties](https://github.com/jonnyhuck/globe-gore-evauation/blob/main/out/12_Gores.png)

And the mathematical 'fit' of the globe to the sphere (shown here for a 12-gore globe):

![fit properties](https://github.com/jonnyhuck/globe-gore-evauation/blob/main/out/12_gore_width.png)

For most purposes, this analysis demonstrates that the [**Cassini Projection**](https://en.wikipedia.org/wiki/Cassini_projection) is best suited to globe production. 