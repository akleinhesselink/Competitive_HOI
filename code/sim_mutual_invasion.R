library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(parallel)

rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

load('data/parms.rda')
load('data/experiments.rda')
load('data/eqs.rda')
load('data/population_sims.rda')

t <- 500 

eqs <- split(eqs, 1:nrow(eqs))
invasion_experiments <- lapply( eqs, function(x){ x[ x == 0 ] <- 0.01 ; x } )
invasion_experiments <- lapply(invasion_experiments, as.numeric)

invasion_results <- lapply( invasion_experiments, run_multi_gen, t = t, parms = parms, tol = 1e-5)

ends <- lapply( invasion_results, function(x){ x <- x[complete.cases(x), ] ; x[nrow(x), ] }  )

invasion_outcomes <- do.call( rbind, ends)
invasion_experiments <- do.call(rbind, invasion_experiments)

saveRDS(invasion_outcomes, 'data/invasion_outcomes.rda')
saveRDS(invasion_experiments, 'data/invasion_experiments.rda')
