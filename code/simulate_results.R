rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
library(grid)

source('code/model_functions.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200              # length of simulation in days 
R_init <- soil_m <- 200   # initial soil moisture (mm water in upper 500 mm of soil)
r <- c(4.2, 2.6, 2.1)     # max uptake rates mm of water per g of plant per day
K <- c(150, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                 # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1         # proportion live biomass converted to seed mass 
R <- seq(0, 200, length.out = 1000)
parms <- list( r = r, 
               K = K, 
               m = m, 
               q = q, 
               soil_m = soil_m, 
               conversion = conversion, 
               seedling_mass = seedling_mass, 
               R = R, 
               times = times)

save(parms,file = 'data/parms.rda')

# Run response surface experiments --------------------------- # 

state <- c( R_init, 0, 0, 0)

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
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  state[4] <- B_init[i,3]*parms$seedling_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, 200, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  
}

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

save(sim_results, file = 'data/sim_results.rda')


