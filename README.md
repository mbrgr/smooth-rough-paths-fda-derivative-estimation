# smooth rough paths fda derivative estimation

This repository contains all numerical results and code implementations associated with the paper:

# Smooth and Rough Paths in Mean Derivative Estimation for Functional Data
  
  Max Berger, Hajo Holzmann
  [arXiv:2503.24066](https://arxiv.org/abs/2503.24066)

##Overview

The paper investigates the estimation of partial derivatives of the mean function in multivariate functional data analysis. It derives near-optimal convergence rates in the minimax sense under a fixed synchronous design over Hölder smoothness classes. 
The study emphasizes the supremum norm due to its relevance in visualizing estimation errors and constructing uniform confidence bands. A key finding is that the smoothness of the process paths significantly influences the rates of convergence:

  - Smooth Paths: If the process paths exhibit higher-order smoothness than the derivative order, the parametric $\sqrt{n}$ rate is achievable under a sufficiently dense design.

  - Rough Paths: For processes with lower-order smoothness, convergence rates are necessarily slower, and the paper determines a near-optimal rate at which estimation remains feasible.

The paper proposes a method to evaluate the smoothness of sample paths by comparing restricted estimates of the partial derivatives of the covariance kernel. 
Therefore the R-Package [biLocPol](https://github.com/mbrgr/biLocPol) (still in work) can be used, which can be installed via:

``` r
# install.packages("devtools")
devtools::install_github("mbrgr/biLocPol")
# library(biLocPol)
```
 

## Repository Contents

The repository is organized as follows:

  - functions.R: Core functions for local polynomial derivative estimation.
  - Figure_1_illustration.R: Script to reproduce Figure 1, illustrating key concepts.
  - Figure_2_optimal_bw_derivative.R: Script for Figure 2, showcasing the sup error with different bandwidths for derivative estimation.
  - Figure_3_4_derivative_error_decomposition.R: Script for Figures 3 and 4, detailing derivative error decomposition.
  - Figure_14_bandwidth_comp_per_degree_full_interval.R: Script for Figure 14, comparing bandwidths per degree over the full interval (analogue to Figure 2)
  - Figure_15_bandwidth_comp_per_degree_full_interval_small_errors.R: Script for Figure 15, full interval, small error variance. (analogue to Figure 2 and 14)
  - Table_1_smooth_and_rough_paths.R: Script to generate Table 1, illustrating the rates of convergence for smooth and rough processes. 
  - Biomechanics.R: Analysis script for the biomechanics dataset. For the analysis of the covariance the biLocPol packages is needed.
  - Weather_in_Nuremberg.R: Analysis script for the weather dataset from Nuremberg. For the analysis of the covariance the biLocPol packages is needed.
  - data/: Includes the results for the simulations and the weather data from Nuremberg.
  - graphics/: Directory for storing generated figures and plots.
  - smooth-rough-paths-fda-derivative-estimation.Rproj: R project file for easy project setup
