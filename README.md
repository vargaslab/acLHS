# An autocorrelated conditioned Latin hypercube method for data-driven spatio-temporal sampling and predictions

by Van Huong Le<sup>1</sup>, Rodrigo Vargas<sup>1</sup>

<sup>1</sup>Department of Plant and Soil Sciences, University of Delaware, Newark, DE, 19716, USA

Corresponding author affiliation and e-mail:

Rodrigo Vargas

Department of Plant and Soil Sciences, University of Delaware, Newark, DE, 19716, USA

[rvargas\@udel.edu](mailto:rvargas@udel.edu)

## Description

This repository contains the source code to perform different sampling methods:Fixed sampling, conditioned Latin Hypercube Sampling (cLHS), and autocorrelated conditioned Latin Hypercube Sampling (acLHS). The predictions are then based on the samples using a Bernstein copula-based stochastic cosimulation (BCSCS) method. Here are two case studies using data of soil CO$_{2}$ efflux (i.e., the CO$_{2}$ efflux from soils to the atmosphere known as soil respiration) that are relevant for carbon cycle science. The first case study represents data from a time series (1D approach), and the second case represents spatial data (2D approach) across the conterminous United States (CONUS).

## Content

##### RProject_acLHS_1D: The case study represents data from a time series

-   Data: this folder contains the data
-   Functions: this folder contains the useful functions
-   Scripts: this folder contains the scripts
-   Results: this folder contains the results
-   RProject_acLHS_1D.Rproj: This file is the R Project

##### RProject_acLHS_2D: The case study represents data from spatial data

-   Data: this folder contains the data
-   Functions: this folder contains the useful functions
-   Scripts: this folder contains the scripts
-   Results: this folder contains the results
-   RProject_acLHS_2D.Rproj: This file is the R Project

## install

The code has been tested using packages of:

-   R version 4.2.1

-   RStudio 2022.07.1

## How to run the code?

#### RProject_acLHS_1D

Opening the project `RProject_acLHS_1D.Rproj` with Rstudio. Then open all the scripts in the "scripts" folder. The scripts are run in the following order: 0_Getting_Started.R, 1_Exploratory_data_analysis.R, 2_Variogram_analysis.R, 3_Sampling_Design.R, 4_Simulations.R.

-   0_Getting_Started.R: this script is for installing and loading R packages, and also loading functions from the functions folder.
-   1_Exploratory_data_analysis.R: this script is for exploring and calculating univariate statistical properties and dependency relationships between variables.
-   2_Variogram_analysis.R: this script is to explore the temporal or spatial distribution of the variable of interest and calculate its autocorrelation function.
-   3_Sampling_Design.R: this script is for applying sampling methods based on the data.
-   4_Simulations.R: this script is to model the characteristic functions of the variables and perform the simulation.

#### RProject_acLHS_2D

Opening the project `RProject_acLHS_2D.Rproj` with Rstudio. Then open all the scripts in the "scripts" folder. The scripts are run in the following order: 0_Getting_Started.R, 1_Exploratory_data_analysis.R, 2_Variogram_analysis.R, 3_Sampling_Design.R, 4_Simulations.R.

-   0_Getting_Started.R: this script is for installing and loading R packages, and also loading functions from the functions folder.
-   1_Exploratory_data_analysis.R: this script is for exploring and calculating univariate statistical properties and dependency relationships between variables.
-   2_Variogram_analysis.R: this script is to explore the temporal or spatial distribution of the variable of interest and calculate its autocorrelation function.
-   3_Sampling_Design.R: this script is for applying sampling methods based on the data.
-   4_Simulations.R: this script is to model the characteristic functions of the variables and perform the simulation.

## License

MIT License

Copyright (c) 2022 Van Huong Le, Rodrigo Vargas
 
