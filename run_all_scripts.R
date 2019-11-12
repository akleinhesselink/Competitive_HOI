#
# Running this script will run all analyses and generate all figures in the 
# manuscript "Mechanisms underlying higher order interactions: from definitions 
# to ecological processes".  Running this script may take over 15 minutes. 
#
# Before running this script set the working directory to the directory containing
# the README (the directory this script is in). 
# 

rm(list = ls())

library( tidyverse ) 
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
