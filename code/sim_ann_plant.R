library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)
rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# set parameters ------------------------------------- 
alphas <- matrix( c(1, 0.5, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 20), 3, 3, byrow = T)
lambdas <- c(24, 32, 41)
taus <- c(-1, -1, -0.2)
pars <- list( lambdas = lambdas, alphas = alphas , taus = taus) 

# ------------------------------------------------ 

experiments <- expand.grid( N1 = c(0, 2^seq(0, 12, 1)), N2 = c(0, 2^seq(0, 12, 1)), N3 = c(0, 2^seq(0, 12, 1)) )
form <- as.formula('~ -1 + N1 + N2 + N3')

out <- experiments
for( i in 1:nrow(experiments)){ 
  seeds <- experiments[i, ]
  out[i, ] <- ann_plant_mod(seeds, form, unlist(pars))
}

experiments <- data.frame(experiments, out)
names(experiments) <- c(paste0('N', 1:3), paste0('F', 1:3))

experiments <- 
  experiments %>% 
  filter( (N1 == 0 & N2 == 0) | (N3 == 0 & N2 == 0 ) | (N1 == 0 & N3 == 0 ) ) %>% 
  mutate( lambda =  ifelse(N1 == 0 & N2 == 0 & N3 == 0 , T, F)) %>% 
  gather( competitor, density, N1:N3) %>% 
  filter( lambda | density > 0 ) %>% 
  gather( focal, fecundity, F1:F3)  


test <- experiments
test$competitor_label <- paste0( 'competitor\n', test$competitor) 
test$focal_label <- paste0( 'focal\n', str_replace( test$focal, 'F', 'N'))

test %>% 
  ggplot(aes( x = density, y = fecundity, color = competitor_label) ) + 
  geom_point() + 
  geom_line() + 
  facet_grid(focal_label ~ competitor_label )


fit1 <- fit_ann_plant(experiments, focal = 1)
fit2 <- fit_ann_plant(experiments, focal = 2)
fit3 <- fit_ann_plant(experiments, focal = 3)

fit1$par
fit2$par
fit3$par

matrix( c(fit1$par[2:4] , fit2$par[2:4], fit3$par[2:4]), 3, 3, byrow = T)
alphas


