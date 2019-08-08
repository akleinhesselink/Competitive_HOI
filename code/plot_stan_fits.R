rm(list = ls())
library(rstan)
library(tidyverse)
library(stringr)
library(gridExtra)
library(grid)

source('code/plotting_parameters.R')
source('code/simulation_functions.R')

load( 'output/parms.rda')
load('output/stan_dat.rda')
load('output/single_sp_fits.rda')
load('output/HOI_fits.rda')


get_pars <- function(x, pars ){ 
  
  pars[str_detect(names(pars), x)]
}

bh_comp <- function(x, pars){ 
  lambda <- get_pars('lambda', pars)
  alpha  <- get_pars('alpha', pars)
  tau    <- get_pars('tau', pars)

  alpha <- matrix( get_pars('alpha', pars), 3, 3, byrow = T)
  
  lambda/(1 + alpha%*%x)^tau
}

bh2_comp <- function(x, pars){ 
  lambda <- get_pars('lambda', pars)
  alpha <- matrix( get_pars('alpha', pars), 3, 3, byrow = T)
  tau   <- matrix( get_pars('tau', pars), 3, 3, byrow = T)  
  y <- NA
  
  for( i in 1:3){ 
    y[i] <- lambda[i]/(1 + sum( (alpha[i, ]*x)^tau[i, ] ) )
  }
  return(y)
}

ricker_comp <- function(x, pars){ 
  lambda <- get_pars('lambda', pars)
  alpha <- matrix( get_pars('alpha', pars), 3, 3, byrow = T)
  tau   <- matrix( get_pars('tau', pars), 3, 3, byrow = T)  
  
  y <- NA
  for( i in 1:3){ 
    y[i] <- lambda[i]*exp(-sum( (alpha[i, ]*x)^tau[i, ] ) )
  }
  return(y)
}

bh2_HOI_comp <- function(x, h, pars){ 
  lambda <- get_pars('lambda', pars)
  alpha <- matrix( get_pars('alpha', pars), 3, 3, byrow = T)
  tau   <- matrix( get_pars('tau', pars), 3, 3, byrow = T)  
  beta  <- get_pars('beta', pars)
  
  y <- NA
  
  for( i in 1:3){ 
    y[i] <- lambda[i]/(1 + sum( (alpha[i, ]*x)^tau[i, ] ) + beta[i]*sqrt(h) )
  }
  return(y)
}

ricker_HOI_comp <- function(x, h, pars){ 
  lambda <- get_pars('lambda', pars)
  alpha <- matrix( get_pars('alpha', pars), 3, 3, byrow = T)
  tau   <- matrix( get_pars('tau', pars), 3, 3, byrow = T)  
  beta  <- get_pars('beta', pars)
  
  y <- NA
  
  for( i in 1:3){ 
    y[i] <- lambda[i]*exp( - sum( (alpha[i, ]*x)^tau[i, ] ) - beta[i]*sqrt(h) )
  }
  return(y)
}


par_ests <- function( fit, parnames){ 
  summary(fit, parnames)$summary[, 1]
}

parnames <- c('lambda', 'alpha', 'tau')

bh_pars <- par_ests(bh_fit, parnames)
bh2_pars <- par_ests(bh_fit2, parnames)
ricker_pars <- par_ests(ricker_fit, parnames)

x <- stan_dat_1c$x[stan_dat_1c$focal == 1]
focal <- stan_dat_1c$focal
y <- stan_dat_1c$y

bh2 <- as.vector( do.call( rbind, lapply( lapply(x, bh2_comp, bh2_pars), t)))
bh <- as.vector( do.call( rbind, lapply( lapply(x, bh_comp, bh_pars), t)))
rick <- as.vector( do.call( rbind, lapply( lapply(x, ricker_comp, ricker_pars), t)))

test <- data.frame(focal = focal, 
           x = do.call(rbind, x),
           y = y, 
           bh = bh, 
           bh2 = bh2, 
           rick = rick)

plot(test$y, test$bh)
plot(test$y, test$bh2)
plot(test$y, test$rick)

parnames <- c(parnames, 'beta')
bh_HOI_pars <- par_ests(bh_HOI_fit, parnames)
ricker_HOI_pars <- par_ests(ricker_HOI_fit, parnames)

x <- stan_dat_2c$x[ stan_dat_2c$focal == 1 ]
y <- stan_dat_2c$y

xmat <- data.frame( do.call(rbind, x))
h <- rowSums( model.matrix(data = xmat , ~ -1 + X1*X2*X3 )[, 3:5])

focal <- stan_dat_2c$focal
bh2_HOI <- as.vector( do.call( rbind, lapply( mapply(FUN = bh2_HOI_comp, x, h, MoreArgs = list( 'pars' = bh_HOI_pars), SIMPLIFY = F), t )) )

data.frame( do.call( rbind, x ) , h)



lambda <- bh_HOI_pars[str_detect(names(bh_HOI_pars), 'lambda') ]
alpha <- bh_HOI_pars[str_detect(names(bh_HOI_pars), 'alpha')]
beta <- bh_HOI_pars[str_detect(names(bh_HOI_pars), 'beta')]
tau <- bh_HOI_pars[str_detect(names(bh_HOI_pars), 'tau')]

HOI_frame <- data.frame( focal = stan_dat_2c$focal, y = stan_dat_2c$y )
HOI_frame <- cbind( HOI_frame, as.data.frame( do.call( rbind, stan_dat_2c$x)) )

mygrid <- expand.grid(focal = c(1:3), 
                      V1 =seq(0,8, by = 0.5), 
                      V2 = seq(0,8, by = 0.5), 
                      V3 = seq(0,8, by = 0.5)) %>% 
          filter( (V1 > 0) + (V2 > 0) + (V3 > 0) < 3) %>%  # drop rows with three species
  left_join(HOI_frame)

test <- 
  HOI_frame %>% 
  ungroup() %>% 
  mutate( trial = row_number()) %>% 
  rename( 'species' = focal) %>% 
  rename( 'B1' = V1, 'B2' = V2, 'B3' = V3) %>% 
  mutate( species_lab = factor(species, labels = c('A) Early', 'B) Mid', 'C) Late')))



test <- 
  test %>% 
  gather( mod_type, y_hat, m4:m5 ) %>% 
  group_by(species, mod_type) %>% 
  mutate( lambda = max(y)) %>% 
  mutate( lambda_plot = ifelse(y == max(y), T, F)) 

theme2 <- 
  journal_theme + 
  theme(strip.text = element_text(hjust = 0.1), 
        legend.title = element_blank(), 
        legend.background = element_blank(), 
        legend.justification = c(1, 1), 
        legend.position = c(1,0.8))

p1 <- 
  test %>% 
  filter( B2 < 15, B3 < 15) %>% 
  filter( species == '1' , B1 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B2, y = y, color = factor(B3), shape = factor(B3) )) + 
  geom_point(size = 3)  + 
  geom_line(aes( y = y_hat, linetype = mod_type)) +
  scale_x_continuous(breaks = c(0,4,8)) + 
  scale_color_manual(values = c('black', 'blue', 'red')) + 
  scale_shape_manual(values = c(1, 0, 2)) + 
  xlab('Mid Competitor\nDensity') + 
  ylab( 'Per Capita Seed Production') + 
  facet_wrap(~ species_lab) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1)) +
  theme2 + 
  annotate(geom = 'text', 
           x = 8, 
           y = 55, 
           label = 'Late Competitor\nDensity', 
           size = 5, 
           hjust = 1)
p1

p2 <- 
  test %>% 
  filter( B1 < 15, B3 < 15) %>% 
  filter( species == '2' , B2 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B3), shape = factor(B3)) ) + 
  geom_point(size = 3) + 
  geom_line(aes( y = y_hat, linetype = mod_type)) + 
  scale_x_continuous(breaks = c(0,4,8)) + 
  scale_color_manual(values = c('black', 'blue', 'red')) + 
  scale_shape_manual(values = c(1, 0, 2)) + 
  xlab('Early Competitor\nDensity') + 
  ylab( 'Per Capita Seed Production') + 
  facet_wrap(~ species_lab) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1)) +
  theme2 + 
  theme(axis.title.y = element_blank()) + 
  annotate(geom = 'text', 
           x = 8, 
           y = 76, 
           label = 'Late Competitor\nDensity', 
           size = 5, 
           hjust = 1)

p2

p3 <- 
  test %>% 
  filter( B1 < 15, B2 < 15) %>% 
  filter( species == '3' , B3 == 0 ) %>% 
  filter( B2 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B2), shape = factor(B2) ) ) + 
  geom_point(size = 3) + 
  geom_line(aes( y = y_hat, linetype = mod_type)) + 
  scale_x_continuous(breaks = c(0,4,8)) + 
  scale_color_manual(values = c('black', 'orange', 'red')) + 
  scale_shape_manual(values = c(1, 0, 2)) + 
  xlab('Early Competitor\nDensity') + 
  ylab( 'Per Capita Seed Production') + 
  facet_wrap(~ species_lab) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1)) +
  theme2 + 
  theme(axis.title.y = element_blank()) + 
  annotate(geom = 'text', 
           x = 8, 
           y = 100, 
           label = 'Mid Competitor\nDensity', 
           size = 5, 
           hjust = 1)

p3

n_comp2 <- grid.arrange(p1, p2, p3, 
                        nrow = 1, 
                        widths = c(0.32, 0.3, 0.3), 
                        top = textGrob('Focal Species', gp = gpar(fontsize = 16)))

