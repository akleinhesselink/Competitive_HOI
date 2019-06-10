rm(list = ls())

library(rstan)
library(tidyverse)
library(stringr)
library(gridExtra)
library(grid)

rstan_options(auto_write = T)

source('code/plotting_parameters.R')
source('code/simulation_functions.R')

load( 'output/sim_results.rda')
load( 'output/parms.rda')

# convert biomass results into final seed production per plant 
sim_results <- 
  sim_results %>% 
  mutate( Y1 = parms$conversion*X2/parms$seed_mass/B1, 
          Y2 = parms$conversion*X3/parms$seed_mass/B2, 
          Y3 = parms$conversion*X4/parms$seed_mass/B3) %>% 
  select( - X1)

sim_results <- 
  sim_results %>% 
  gather( species, y, Y1, Y2, Y3)  %>% 
  mutate( y = y + rnorm(1, 0, 0.001)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( B3 = ifelse( str_extract(species, '\\d+') == 3 & B3 > 0, B3 - 1, B3)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) ) %>% 
  filter( n_comp < 3)


# add HOI column to data ---------------- # 

sim_results <- sim_results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) )

# Filter out NA ------------------------- #  

sim_results <- sim_results %>% 
  filter( !is.na(y))

all_dat <- 
  sim_results %>% 
  distinct() %>% 
  rename( 'focal' = species )  %>% 
  group_by( row_number()) %>% 
  mutate( x = list( c(B1, B2, B3))) %>% 
  ungroup %>% 
  mutate( focal = as.numeric(factor(focal)))  %>% 
  select( focal, y, n_comp, x)

make_stan_dat <- 
  function(df ,  subset_cond ) { 
  
  df <- df %>% 
    filter( eval( parse( text = subset_cond ) ))
  
  list( 
    focal = df$focal, 
    y = df$y,
    x = df$x,
    N = nrow(df), 
    S = length(unique( df$focal))
    )
}

stan_dat_1c <- make_stan_dat(all_dat, 'n_comp < 2' )
stan_dat_2c <- make_stan_dat(all_dat, 'n_comp < 3' )

# make smoother prediction grid for plotting curves
smooth_grid <- expand.grid( focal = 1:3, x1 = seq(0, 16, by = 0.1), x2 = seq(0, 16, by = 0.1), x3 = seq(0, 16, by = 0.1))
smooth_grid$x <- as.matrix( smooth_grid[, c('x1', 'x2', 'x3')] )
smooth_grid_1c <- smooth_grid[ rowSums(smooth_grid$x > 0 ) < 2, ]
smooth_grid_2c <- smooth_grid[ rowSums(smooth_grid$x > 0 ) < 3, ]


save(file = 'output/stan_dat.rda', list = c('stan_dat_1c', 'stan_dat_2c', 'smooth_grid_1c', 'smooth_grid_2c'))
