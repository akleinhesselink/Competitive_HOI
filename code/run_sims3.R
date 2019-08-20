rm(list = ls())

library(deSolve)
library(tidyverse)
source('code/simulation_functions3.R')
load( 'output/parms3.rda')

max_d <- 36
min_d <- 0 
dgrad <- sort( c(2,3,(seq(sqrt(min_d), sqrt(max_d), 1))^2) )

B_init <- expand.grid( 
  B1 = dgrad, 
  B2 = dgrad,
  B3 = dgrad)

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 ) %>%  # filter out three species cases 
  filter( B1 + B2 + B3 < 50) # filter out cases where total density is high  

state <- c( R = parms$R0, B1 = 0, B2 = 0, B3 = 0)

out <- list()

for( i in 1:nrow(B_init)){
  parms$n <- as.numeric( B_init[i, 1:3] )
  parms$n[parms$n == 0 ] <- 1
  
  state[2] <- B_init[i,1]*parms$seed_mass
  state[3] <- B_init[i,2]*parms$seed_mass
  state[4] <- B_init[i,3]*parms$seed_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, parms$U, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root3, 
                  event = list(func = event, root = T))
  
}

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

sim_results <- 
  sim_results %>% 
  rename( "X1" = R, "X2" = B1.1, "X3" = B2.1, "X4" = B3.1)

save(sim_results, file = 'output/sim_results3.rda')

