rm(list = ls())

library(deSolve)
library(tidyverse)
source('code/simulation_functions3.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
U <- 200                    # Length of simulation in days 
R0 <- 400                   # Initial resource concentration 
Vmax <- 0.12                # Max R uptake rates per unit surface area 
b0 <- 0.3                   # Base of resource uptake function 
K <- 200                    # Resource half-saturation constants
nu <- 0.66                  # Exponent for scaling of tissue mass to surface area 
d <- c(0.001, 0.005, 0.01)  # "Density", i.e. investment in resource acquisition tissues per unit mass  
m <- c(0.5, 0.13, 0.07)     # Biomass loss rate per unit mass 
q <- 0.2                    # Resource use efficiency 
seed_mass <- c(0.005)       # Seed/seedling mass 
conversion <- 0.1           # Proportion live biomass converted to seed mass 

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
               U = U)

save(parms, file = 'output/parms3.rda')

# Run response surface experiments --------------------------- # 

state <- c( R = R0, B1 = 0, B2 = 0, B3 = 0)

B_scen <- expand.grid( B1 = c(0, 1), B2 = c(0,1), B3 = c(0,1))
B_scen <- B_scen[c(2, 3, 5, 8), ]

run_basic_scenarios <- function(scenarios, parms){ 
  out <- list() 
  
  for( i in 1:nrow(scenarios)){ 
    parms$n <- as.numeric( scenarios[i, 1:3] )
    parms$n[parms$n == 0 ] <- 1

    state[2] <- scenarios[i,1]*parms$seed_mass
    state[3] <- scenarios[i,2]*parms$seed_mass
    state[4] <- scenarios[i,3]*parms$seed_mass 

    out[[i]] <- ode(state, 
                times = seq(1, parms$U, by = 0.1), 
                func = grow, 
                parms = parms, 
                rootfun = root3, 
                event = list(func = event, root = T))
  }
  return(out)
}

basics <- run_basic_scenarios(B_scen, parms)

labels <- c('early', 'mid', 'late')
cols <- c('black', 'red', 'blue')

show_results <- function( out, labels, cols ){ 
  
  par(mfrow = c(length(out),1), mar = c(2,5,1,1))
  for(i in 1:(length(out)-1) ){ 
    plot(out[[i]][, 1], out[[i]][, i+2], type = 'l', ylab = labels[i], 
         col = cols[i], lty = 1)
  }
  par(mar = c(2,5,1,1), xpd = T)

  matplot(out[[length(out)]][, 1], 
          out[[length(out)]][, c((length(out) - 1):(length(out) + 1))], 
          type = 'l', 
          ylab = 'Biomass', 
          xlab = 'Days', col = cols, lty = 1)
  y_pos <- 0.8*max( apply( out[[4]][, 3:(length(out) + 1)], 2, max) )
  legend(x = 160, y = y_pos, legend = labels, col = cols, lty = 1, cex = 1.2)

}

show_results(basics, labels, cols)

par(mfrow = c(1,1), mar = c(4,4,3,3))
plot( d, m, type = 'b', col = cols, ylab = 'Loss rate (m)', xlab = 'Tissue density (d)')
text( d[1], m[1], labels[1], adj = c(-0.5, 0.5))
text( d[2], m[2], labels[2], adj = c(-0.5, -1))
text( d[3], m[3], labels[3], adj = c(0.5, -1))

lm(data =  data.frame(d = d, m = m) , d ~ m )
