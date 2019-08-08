rm(list = ls())

library(rstan)
library(tidyverse)
library(stringr)
library(gridExtra)
library(grid)

rstan_options(auto_write = T)

source('code/plotting_parameters.R')
source('code/simulation_functions.R')

load('output/parms.rda')
load('output/stan_dat.rda')

N <- stan_dat_2c$N
combos <- t (combn(3:1, 2) )
x <- do.call( rbind, stan_dat_2c$x )
H <- x 
H[] <- NA

for( i in 1:N) { 
  for( j in 1:nrow(combos)){
    H[i, j] <- x[i, combos[j, 1]]*x[i, combos[j, 2]]
  }
}

stan_dat_2c$h <- rowSums(H) # collapse to a single term giving interspecific HOI 




bh_HOI_fit <- stan(file = 'code/bh_comp_HOI.stan', 
            data = stan_dat_2c, 
            chains = 4, 
            cores = 4)

ricker_HOI_fit <- stan(file = 'code/ricker_comp_HOI.stan', 
                  data = stan_dat_2c, 
                  chains = 4, 
                  cores = 4, 
                  control = list('adapt_delta' = 0.9))

stan_dat_2c_with_h <- stan_dat_2c

save(bh_HOI_fit, ricker_HOI_fit, stan_dat_2c_with_h, file = 'output/HOI_fits.rda')


