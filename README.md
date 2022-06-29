# VStuning_project
This repository is for the analysis of VS neuron recordings developped in the [Bagnall Lab](https://sites.wustl.edu/bagnall/). All codes are written in Matlab

Please see specifics in our paper: [Central vestibular tuning arises from patterned convergence of otolith afferents.](https://www.biorxiv.org/content/10.1101/2020.02.14.948356v1)


1. The event detection algorithm, used a derivative-based method to detect EPSCs and EPSPs
The algorithm was adapted from a previous method described in:
[Bagnall MW, McLean DL. Modular organization of axial microcircuits in zebrafish. Science. 2014 Jan 10;343(6167):197-200.](https://science.sciencemag.org/content/343/6167/197.full)

2. The deconvolution algorithm used a sparse deconvolution method based on [FISTA](https://github.com/tiepvupsu/FISTA), as well as [ISO-SPILT](https://github.com/flatironinstitute/isosplit5).
