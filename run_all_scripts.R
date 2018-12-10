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

source('code/run_simulations.R')
source('code/make_figure_two.R')
source('code/fit_models.R')
source('code/make_figure_six.R')
source('code/vary_tradeoff.R')

