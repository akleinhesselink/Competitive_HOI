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

traceplot(out, 'tau')

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
tau_fit <- matrix( tau_fit, 3,3, byrow = T)

alpha_fit <- par_ests[grep('alpha', names(par_ests))]
alpha_fit <- matrix( alpha_fit, 3, 3 , byrow = T)


# Predict Y on smooth grid  ----------------------- # 
S <- length(unique( smooth_grid_1c$focal))
N <- nrow(smooth_grid_1c)
test_y_hat <- NA
for( i in 1:N) {
  C <- NA
  for(j in 1:S){ 
    C[j] <- (alpha_fit[smooth_grid_1c$focal[i], j]*smooth_grid_1c$x[i, j])^tau_fit[smooth_grid_1c$focal[i], j]
  }
  test_y_hat[i] <- lambda_fit[smooth_grid_1c$focal[i]]/(1 + sum(C))
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


# Predict Y on 2-d smooth grid (multiple competitors) # 

S <- length(unique( smooth_grid_2c$focal))
N <- nrow(smooth_grid_2c)
test_y_hat <- NA
for( i in 1:N) {
  C <- NA
  for(j in 1:S){ 
    C[j] <- (alpha_fit[smooth_grid_2c$focal[i], j]*smooth_grid_2c$x[i, j])^tau_fit[smooth_grid_2c$focal[i], j]
  }
  test_y_hat[i] <- lambda_fit[smooth_grid_2c$focal[i]]/(1 + sum(C))
}

smooth_grid_2c$y_hat <- test_y_hat

# species 1 -------------------------- # 

pred_surface1 <- 
  smooth_grid_2c %>% 
  select( - x) %>% 
  arrange(focal, x1, x2, x3) %>% 
  filter( focal == 1, x1 == 0 )  %>% 
  select( x2, x3, y_hat) %>% 
  spread(x3, y_hat) %>% 
  as.matrix()

pred_surface1 <- pred_surface1[, -1]

# species 2 ---------------------------- # 

pred_surface2 <- 
  smooth_grid_2c %>% 
  select( - x) %>% 
  arrange(focal, x1, x2, x3) %>% 
  filter( focal == 2, x2 == 0 )  %>% 
  select( x1, x3, y_hat) %>% 
  spread(x3, y_hat) %>% 
  as.matrix()

pred_surface2 <- pred_surface2[, -1]

# species 3 ---------------------------- # 

pred_surface3 <- 
  smooth_grid_2c %>% 
  select( - x) %>% 
  arrange(focal, x1, x2, x3) %>% 
  filter( focal == 3, x3 == 0 )  %>% 
  select( x1, x2, y_hat) %>% 
  spread(x2, y_hat) %>% 
  as.matrix()

pred_surface3 <- pred_surface3[, -1]

save(file = 'output/pred_surfaces.rda', list = c('pred_surface1', 'pred_surface2', 'pred_surface3'))
