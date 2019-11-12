rm(list = ls())

library(deSolve)
library(tidyverse)
source('code/simulation_functions2.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
U <- 200                  # length of simulation in days 
R0 <- 200                 # initial resource concentration 
r <- c(4.0, 3.9, 3.8)     # max uptake rates of resource per unit mass of plant per day
K <- c(50, 50, 50)        # resource half-saturation constants
m <- c(0.09, 0.09, 0.09)  # tissue respiration and loss rate units of mass per unit plant mass per day 
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

eventdat <- data.frame(var = c("R", "B1", "B2", "B3"),
                       time = c(rep(50, 4), rep(60, 4), rep(70, 4)),
                       value = c(c(1,0,1,1), c(1,1,0,1), c(1,1,1,0)),
                       method = c("mult", "mult", "mult"))

parms$eventdat <- eventdat

save(parms, file = 'output/parms2.rda')

# Run response surface experiments --------------------------- # 

state <- c( R = R0, B1 = 0, B2 = 0, B3 = 0)

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
                  event = list(data = eventdat))
  
}

get_last <- function(x) { 
  y <- tail(x[x > 0], 1)  
  ifelse( length(y) == 0 , 0, y)  
}

sim_results2 <- do.call( rbind, lapply( out, function(x) apply( x, 2, get_last )))
sim_results2 <- data.frame(B_init, sim_results2 )

sim_results2 <- 
  sim_results2 %>% 
  rename( "X1" = R, "X2" = B1.1, "X3" = B2.1, "X4" = B3.1)


save(sim_results2, file = 'output/sim_results2.rda')

