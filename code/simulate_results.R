rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
library(grid)

source('code/model_functions.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
U <- 200                  # length of simulation in days 
R0 <- 200                 # initial resource concentration 
r <- c(4.2, 2.6, 2.1)     # max uptake rates of resource per unit mass of plant per day
K <- c(150, 30, 0.5)      # resource half-saturation constants
m <- 0.09                 # tissue respiration and loss rate units of mass per unit plant mass per day 
q <- 0.07                 # resource use efficiency 
seed_mass <- c(0.005)     # seed/seedling mass 
conversion <- 0.1         # proportion live biomass converted to seed mass 

parms <- list( r = r, 
               K = K, 
               m = m, 
               q = q, 
               R0 = R0, 
               conversion = conversion, 
               seed_mass = seed_mass, 
               #R = R, 
               U = U)

save(parms,file = 'data/parms.rda')

# Run response surface experiments --------------------------- # 

state <- c( R0, 0, 0, 0)

B_init <- expand.grid( 
  B1 = c(0, seq(1, 8, by = 1), 16), 
  B2 = c(0, seq(1, 8, by = 1), 16), 
  B3 = c(0, seq(1, 8, by = 1), 16))

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 )  # filter out three species cases 

out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]*parms$seed_mass
  state[3] <- B_init[i,2]*parms$seed_mass
  state[4] <- B_init[i,3]*parms$seed_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, parms$U, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  
}

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

save(sim_results, file = 'data/sim_results.rda')


