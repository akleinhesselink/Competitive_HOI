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

bh_fit <- stan(file = 'code/bh_comp.stan', 
             data = stan_dat_1c, 
             chains = 4, 
             cores = 4)

bh_fit2 <- stan(file = 'code/bh_comp2.stan', 
               data = stan_dat_1c, 
               chains = 4, 
               cores = 4)

ricker_fit <- stan(file = 'code/ricker_comp.stan', 
                   data = stan_dat_1c, 
                   chains = 4, 
                   cores = 4)


save(bh_fit, bh_fit2, ricker_fit, file = 'output/single_sp_fits.rda')

par_ests <- summary(out, c('lambda', 'tau', 'alpha'))$summary[, '50%']
y_hat <- summary(out, 'y_hat')
y_hat <- data.frame( y_hat$summary) %>% select( mean, X50.)
test <- data.frame( focal = stan_dat_1c$focal, y = stan_dat_1c$y )
test <- cbind( test, as.data.frame( do.call( rbind, stan_dat_1c$x)) )
test$y_hat <- y_hat$mean

test_long <- 
  test %>% 
  gather( comp, density, V1:V3) %>%
  mutate( comp = as.numeric(factor(comp))) %>% 
  group_by( comp) %>% 
  filter( row_number() == 1 | density > 0)
  
test_long %>% 
  ggplot( aes( x = density, y = y, color = factor( comp ) )) + 
  geom_point() + 
  geom_line(aes( y = y_hat)) + 
  facet_wrap(~focal)

lambda_fit <- par_ests[ grep('lambda', names( par_ests )  ) ] 
tau_fit <- par_ests[grep('tau', names(par_ests))]
alpha_fit <- par_ests[grep('alpha', names(par_ests))]
alpha_fit <- matrix( alpha_fit, 3, 3 , byrow = T)

# Predict y on smooth grid  ----------------------- # 
test_y_hat <- NA
for( i in 1:nrow(smooth_grid_1c)) { 
  test_y_hat[i] <- lambda_fit[smooth_grid_1c$focal[i]]/( 1 + alpha_fit[smooth_grid_1c$focal[i], ]%*%smooth_grid_1c$x[i, ])^tau_fit[smooth_grid_1c$focal[i]]
}
smooth_grid_1c$y_hat <- test_y_hat

smooth_grid_1c <- 
  smooth_grid_1c %>% 
  select( - x) %>% 
  gather( comp, density, x1:x3) %>% 
  mutate( comp = as.numeric(factor(comp)))  %>% 
  arrange(focal, comp, density )  %>% 
  filter( row_number() == 1 | density > 0)  

# ------------------------------------------------ # 

smooth_grid_1c %>% 
  ggplot( aes( x = density, y = y_hat, color = factor(comp) )) + 
  geom_line() + 
  facet_wrap(~focal) + 
  geom_point( data = test_long, aes( y = y)) + 
  scale_color_manual( values = my_colors) + 
  theme_bw()
