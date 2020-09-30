# lmcrCoda-Reproducibility

## About this repository

This repository provides a collection of [R](https://CRAN.R-project.org/) 
scripts to reproduce all examples, simulations and figures from the manuscript:  

N. Stefelova, A. Alfons, J. Palarea-Albaladejo, P. Filzmoser and K. Hron.
Robust regression with compositional covariates including cellwise outliers.
Under review.

Please note that this repository is rather large because it also contains R 
data files with all simulation results.  This way, if you only want to quickly
reproduce the figures with simulation results, you do not actually need to run 
the simulations first.


## Reproduce results

The easiest way to reproduce the results is to clone this repository with 
[RStudio](https://rstudio.com/products/rstudio/download/).  Running the 
scripts within the resulting RStudio project ensures that there are no issues 
with file paths for storing or reading results, or for producing files 
containing plots.  In addition, the RStudio project uses 
[packrat](https://rstudio.github.io/packrat/) to make sure that the correct 
versions of all required R packages are used.  After opening the RStudio 
project for the fist time, please type `packrat::restore()` on the R command 
line to retrieve the correct versions of all required packages.


## R functions included in this repository

The subfolder `code` contains the R script `lmcrCoda.R` with our implementation 
of the cellwise and rowwise robust estimator of linear regression with 
compositional covariates, as well as the R script `shootingS.R` with the 
implementation of the shooting S-estimator (Öllerer, Alfons & Croux, 
2016).  Those functions are only included here for the purpose of reproducing 
the results from the manuscript.

For any other purposes, the most up-to-date versions of the R functions should 
always be obtained from <https://github.com/aalfons/lmcrCoda> and 
<https://github.com/aalfons/shootingS>, respectively.
