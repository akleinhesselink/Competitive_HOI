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

stan_dat <- make_stan_dat(all_dat, 'n_comp < 2' )

out <- stan(file = 'code/bh_comp.stan', 
             data = stan_dat, 
             chains = 4, 
             cores = 4)

traceplot(out, 'tau')

par_ests <- summary(out, c('lambda', 'tau', 'alpha'))$summary[, '50%']

y_hat <- summary(out, 'y_hat')

y_hat <- data.frame( y_hat$summary) %>% select( mean, X50.)

test <- data.frame( focal = stan_dat$focal, y = stan_dat$y )

test <- cbind( test, as.data.frame( do.call( rbind, stan_dat$x)) )

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


smooth_grid <- expand.grid( focal = 1:3, x1 = seq(0, 16, by = 0.1), x2 = seq(0, 16, by = 0.1), x3 = seq(0, 16, by = 0.1))
smooth_grid$x <- as.matrix( smooth_grid[, c('x1', 'x2', 'x3')] )

smooth_grid <- smooth_grid[ rowSums(smooth_grid$x > 0 ) < 2, ]

test_y_hat <- NA
for( i in 1:nrow(smooth_grid)) { 
  test_y_hat[i] <- lambda_fit[smooth_grid$focal[i]]/( 1 + alpha_fit[smooth_grid$focal[i], ]%*%smooth_grid$x[i, ])^tau_fit[smooth_grid$focal[i]]
}

smooth_grid$y_hat <- test_y_hat

smooth_grid <- 
  smooth_grid %>% 
  select( - x) %>% 
  gather( comp, density, x1:x3) %>% 
  mutate( comp = as.numeric(factor(comp)))  %>% 
  arrange(focal, comp, density )  %>% 
  filter( row_number() == 1 | density > 0)

smooth_grid %>% 
  ggplot( aes( x = density, y = y_hat, color = factor(comp) )) + 
  geom_line() + 
  facet_wrap(~focal) + 
  geom_point( data = test_long, aes( y = y)) + 
  scale_color_manual( values = my_colors) + 
  theme_bw()
