22 September 2006


This folder contains 3 files. These files contain a number of R
functions that I have written for ABC and MCMC analysis. Most of the
papers I've written have used some of these functions (or earlier
versions of them). The functions are commented, and (with a bit of
effort) should be straightforward to use.

loc2plot_d.r - contains a number of functions that were originally
used for MCMC analysis (e.g. Beaumont, Genetics, 1999; papers with
Jay Storz; Beaumont, Genetics, 2003) but can be used for ABC
analysis (as in Beaumont et al, Genetics, 2002). They plot bivariate
densities, and get univariate and bivariate summaries of posterior
distributions. These use the density estimation package locfit
(downloadable from CRAN).


make_pd2.r - This contains the function that will do
regression-based ABC, as described in Beaumont et al (Genetics,
2002). It has evolved slightly from the original script used for
that paper, to allow for parameter values to be transformed, and
also to allow for any number of summary statistics (the original
version needed the function to be modified to accommodate the
required number of summary statistics). These last two changes were
based on modifications to the original script made by Shola Ajayi.

calmod.r - this is a script for doing regression-based model choice
as described in:

Beaumont, M.A. (2006). Joint determination of topology, divergence
time, and immigration in population trees. In Simulation, Genetics,
and Human Prehistory, eds. S. Matsumura, P. Forster, & C. Renfrew.
(McDonald Institute Monographs.) Cambridge: McDonald Institute for
Archaeological Research.

calmod uses the VGAM package. You need to get it from the website
(http://www.stat.auckland.ac.nz/~yee).
