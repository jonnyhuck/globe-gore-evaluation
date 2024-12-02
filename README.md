# Evaluating Map Projections for Globe Making

This repository is intended to facilitate the review of the manuscript _Evaluating Map Projections for Globe Making_.

The code for distortion calculations is in `distortions.py` and the actual process of calculating the results (using this file) is undertaken by `evaluate_gores.py`.

The figure illustrating the input gores is generated using `maps_gores.py`.

My implementation of the *Ginzburg & Salmanova (1964) Globe Gore Projection* is in `globe_gore_projection.py`. Note that this is a modification of the equations given in *Bugayevsky & Snyder 1995*, which appears to be incorrect.

The results tables for all of the evaluated combinations of gore number and sphere diameter are given in `results.txt`, the associated figures are in `out/`.

For most purposes, this analysis demonstrates that the [**Cassini Projection**](https://en.wikipedia.org/wiki/Cassini_projection) is best suited to globe production. 

The included data (used for drawing the figure) are from [Natural Earth](https://www.naturalearthdata.com/).