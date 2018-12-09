rm(list = ls())

library( tidyverse ) 
library( stringr )
library( deSolve )
library( grid ) 
library( gridExtra )

if(!dir.exists(paths = 'figures/')){ 
  dir.create('figures/')
}
if(!dir.exists(paths = 'output/')){
  dir.create('output/')
}

source('code/make_figure_two.R')
source('code/run_simulations.R')
source('code/fit_models.R')
source('code/make_mechanism_figure.R')
source('code/vary_tradeoff.R')

