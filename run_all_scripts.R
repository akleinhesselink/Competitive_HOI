#
# Running this script will run all analyses and generate all figures in the 
# manuscript "Detecting and interpreting higher order interactions in 
# ecological communities".  Running this script may take over 15 minutes. 
#
# Either open the R project file in R studio or set the working directory to 
# the directory containing the README (the directory this script is in). 
# 
# Use 'renv::restore()' to install the required package versions if needed. 

renv::restore() # install required packages 

rm(list = ls())

library( tidyr )
library( dplyr )
library( ggplot2 )
library( stringr )
library( deSolve )
library( gridExtra )
library( grid ) 

if(!dir.exists(paths = 'figures/')){ 
  dir.create('figures/')
}
if(!dir.exists(paths = 'output/')){
  dir.create('output/')
}

# Code takes several minutes to run 
source('code/setup_parms.R')
source('code/run_simulations.R')
source('code/fit_models.R')
source('code/plot_fits.R')
source('code/plot_parameters.R')

# For Appendix 
source('code/vary_tradeoff.R') # code takes several minutes to run
source('code/plot_trade_off_results.R')
