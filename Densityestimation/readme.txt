DENSITYESTIMATION.m computes the empirical means (sample means) of Y_i and Y_i Y_j (i~j), that is, the densities of coloured vertices and of adjacent pairs of coloured vertices as a function of lambda (seeding rate) and mu (percolation parameter) when number of colours is 1.
For each parameter vector (lambda, mu), we create a total_nrows x total_ncols lattice, and cut out the middle nrows x ncols sublattice to compute the densities.
One should specify a vector of lambda and another vector of mu values in DENSITYESTIMATION.m, and all combinations will be computed.

Note that this program uses earlier versions of some of the functions that it shares with the main MSM parameter estimation program, e.g. inputs and outputs are different or they come in different orders. Therefore the functions with shared names are not interchangeable between parameter estimation and density estimation.
This implementation could be made faster by integrating CREATESYNTHDATA_RANDOM_M2.m and SIMULATION_M2.m from the parameter estimation program.

Background:
Felix Beck, Bence Mélykúti.
Parameter estimation in a subcritical percolation model with colouring
arXiv:1604.08908
especially Appendix of version 3, arXiv:1604.08908v3
https://arxiv.org/abs/1604.08908

Authors:
Felix Beck, Bence Melykuti (University of Freiburg, Germany)
2015-2017


Structure:

densityestimation
	createsynthdata_determ
	createsynthdata_random
		simulation1/simulation2
			contamination

data_visualise
	data_surface
	data_contour

data_visualise_publication
	data_surface
	data_contour

(data_visualise and data_visualise_publication contain two of their own functions each, data_surface and data_contour)
