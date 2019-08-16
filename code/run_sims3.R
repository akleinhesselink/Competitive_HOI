rm(list = ls())

library(deSolve)
library(tidyverse)
source('code/simulation_functions3.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
U <- 200                    # Length of simulation in days 
R0 <- 200                   # Initial resource concentration 
Vmax <- 0.2                 # Max R uptake rates per unit surface area 
b0 <- 0.3                   # Base of resource uptake function 
K <- 110                    # Resource half-saturation constants
nu <- 0.7                   # Exponent for scaling of tissue mass to surface area 
d <- c(0.0006, 0.005, 0.01) # "Density", i.e. investment in resource acquisition tissues per unit mass  
m <- c(0.7, 0.14, 0.07)     # Biomass loss rate per unit mass 
q <- 0.1                    # Resource use efficiency 
seed_mass <- c(0.005)       # Seed/seedling mass 
conversion <- 0.1           # Proportion live biomass converted to seed mass 

par(mfrow = c(1,1), mar = c(4,4,3,3))
plot( d, m, type = 'b')

parms <- list( Vmax = Vmax, 
               K = K, 
               m = m, 
               q = q, 
               b0 = b0, 
               nu = nu,
               R0 = R0,
               d = d,
               conversion = conversion, 
               seed_mass = seed_mass, 
               #R = R, 
               U = U)

save(parms, file = 'output/parms3.rda')

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
mynums <- c(1, 10, 100, 111)

for( i in 1:nrow(B_init)){
  # for( j in 1:4){
  # i <- mynums[j]
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

# par(mfrow = c(4,1), mar = c(2,3,1,1))
# plot(out[[1]][, 1], out[[1]][, 3], type = 'l')
# plot(out[[2]][, 1], out[[2]][, 4], type = 'l')
# plot(out[[3]][, 1], out[[3]][, 5], type = 'l')
# matplot(out[[4]][, 1], out[[4]][, c(3:5)], type = 'l')

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

sim_results <- 
  sim_results %>% 
  rename( "X1" = R, "X2" = B1.1, "X3" = B2.1, "X4" = B3.1)

save(sim_results, file = 'output/sim_results3.rda')

