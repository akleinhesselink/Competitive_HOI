rm(list = ls())

library(rstan)
library(tidyverse)
library(stringr)
library(gridExtra)
library(grid)

rstan_options(auto_write = T)

source('code/plotting_parameters.R')
source('code/simulation_functions.R')

load( 'output/parms.rda')

load('output/stan_dat.rda')

out <- stan(file = 'code/bh_comp2.stan', 
            data = stan_dat_1c, 
            chains = 4, 
            cores = 4)

y_bh <- summary(out, 'y_hat')
y_bh <- data.frame( y_bh$summary) %>% select( mean)
test <- data.frame( focal = stan_dat_1c$focal, y = stan_dat_1c$y )
test <- cbind( test, as.data.frame( do.call( rbind, stan_dat_1c$x)) )
test$y_bh <- y_bh$mean

test_long <- 
  test %>% 
  gather( comp, density, V1:V3) %>%
  mutate( comp = as.numeric(factor(comp))) %>% 
  group_by( focal,  comp) %>% 
  filter( row_number() == 1 | density > 0)

test_long %>% 
  ggplot( aes( x = density, y = y, color = factor( comp ) )) + 
  geom_point() + 
  geom_line(aes( y = y_bh)) + 
  facet_wrap(~focal)

# fit ricker model 

out <- stan(file = 'code/ricker_comp.stan', 
            data = stan_dat_1c, 
            chains = 4, 
            cores = 4)

y_r <- summary(out, 'y_hat')
y_r <- data.frame( y_r$summary) %>% select( mean)
test$y_r <- y_r$mean

test_long <- 
  test %>% 
  gather(type, y_hat, c('y_r', 'y_bh')) %>% 
  gather( comp, density, V1:V3) %>%
  mutate( comp = as.numeric(factor(comp))) %>% 
  group_by( focal, comp, type) %>% 
  filter( row_number() == 1 | density > 0)

test_long %>% 
  ggplot( aes( x = density, y = y, color = factor( comp ) )) + 
  geom_point() + 
  geom_line(aes( y = y_hat, linetype = type)) + 
  facet_wrap(~focal)

test_long %>% filter( focal == 2, comp == 1)

out <- stan(file = 'code/custom_comp.stan', 
            data = stan_dat_1c, 
            chains = 4, 
            cores = 4)

y_c <- summary(out, 'y_hat')
y_c <- data.frame( y_c$summary) %>% select( mean, X50.)
test$y_c <- y_c$mean

test_long <- 
  test %>% 
  gather(type, y_hat, c('y_r', 'y_bh', 'y_c')) %>% 
  gather( comp, density, V1:V3) %>%
  mutate( comp = as.numeric(factor(comp))) %>% 
  group_by( comp) %>% 
  filter( row_number() == 1 | density > 0)

test_long %>% 
  ggplot( aes( x = density, y = y, color = factor( comp ) )) + 
  geom_point() + 
  geom_line(aes( y = y_hat, linetype = type)) + 
  facet_wrap(~focal)
