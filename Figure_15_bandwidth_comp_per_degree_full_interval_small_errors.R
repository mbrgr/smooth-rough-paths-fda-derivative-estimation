#### description ####
# Figures in the Appendix; analogue to Figure 2 but on the full interval [0,1] and smaller error variance

#### packages ####
library(locpol)
library(tidyverse)
library(biLocPol)
library(parallel)
library(future.apply)

#### source functions ####
source("functions.R")
#### results ####
load("R-Codes/data/bandwidth_comparison_per_degree_full_interval.RData")

##### Bandwidth Comparison: Different degrees  #####
set.seed(264) 
N = 1000
n = 600
p = c(65, 115, 175, 275, 400, 550)
sigma = .1
m = length(p)
H = sapply(p, function(x)(rev(seq(0.3, 4/x, -0.005))))
H


###### Degree = 2 ######
alpha = 2.1
deg = floor(alpha)

cl = makeCluster(detectCores( ) - 1);
plan(multisession);
erg_deg2 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB_derivative(n, p[j], N, H[[j]], deg = deg, rm_boundary_effect = F, 
                             sigma = sigma )) 
})
stopCluster(cl);

erg_deg_2_df = data.frame(sup.err = unlist(erg_deg2), p = 
                            factor(unlist(sapply(1:m, function(j){rep(p[j], length(H[[j]]))}))),
                          h = unlist(H)
) 

###### Figure 15a ######
Figure15a = ggplot(erg_deg_2_df, aes(x = h, y = sup.err, color = p, pch = p)) + 
  geom_point() + labs(subtitle = "n = 600") + lims(y = c(0.5, 20))
Figure15a

ggsave("Grafics/derivative_bandwidth_comparison_quad_full_interval_small_error.png", device = "png", width = 5, height = 3.8, units = "in")

###### Degree = 3 ######
set.seed(264) # same seed as for degree = 2
alpha = 3.1
deg = floor(alpha)

cl = makeCluster(detectCores( ) - 1);
plan(multisession);
erg_deg3 = lapply(1:m, function(j){
  print(j)
  return(h_est_BB_derivative(n, p[j], N, H[[j]], deg = deg, rm_boundary_effect = F, 
                             sigma = sigma )) 
})
stopCluster(cl);

erg_deg_3_df = data.frame(sup.err = unlist(erg_deg3), p = 
                            factor(unlist(sapply(1:m, function(j){rep(p[j], length(H[[j]]))}))),
                          h = unlist(H)
) 

###### Figure 3 ######
# Local cubic estimator bandwidth comparison 
Figure15b = ggplot(erg_deg_3_df, aes(x = h, y = sup.err, color = p, pch = p)) + 
  geom_point() + lims(y = c(0.5, 20)) + labs(subtitle = "n = 600") 
Figure15b

ggsave("grafics/derivative_bandwidth_comparison_cubic_full_interval_small_error.png", device = "png", width = 5, height = 3.8, units = "in")

save.image("data/bandwidth_comparison_per_degree_full_interval_small_error.RData")
