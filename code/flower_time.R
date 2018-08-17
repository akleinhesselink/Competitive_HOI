rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

# graphics themes ------------------------------------------------ # 

journal_theme <- 
  my_theme + 
  theme( axis.title = element_text(size = 12), 
         legend.text = element_text(size = 10), 
         legend.title = element_text(size = 12), 
         strip.text = element_text(size = 12), 
         axis.text = element_text(size = 10))

# Functions ------------------------------------------------------ #

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    R[round(R, 7) == 0 ] <- 0 
    V  <- State[2:4]                             # veg biomass for each species
    A  <- as.numeric( State[5:7] > 0)            # reproduction switch
    
  
    U <- V*(2/3)*f(R, r, K)                          # vector of net resource capture per species
    dR <- - sum(U)                                  # sum of resource uptake from all species 
    
    G <- q*U                                   # vector total growth per species
    
    dV <- NA
    if( A[1] == 0 ) { 
      dV[1] <- G[1] - V[1]*m                              # vector net veg growth 
    
      tlag <- u - 1
      if (tlag < 0)
        dV_lag <- 0
      else 
        dV_lag <- lagderiv(tlag, 2)  
        
      A[1] <- as.numeric(dV_lag > dV[1])
    }
    
    if( A[2] == 0 ) { 
      dV[2] <- G[2] - V[2]*m                              # vector net veg growth 
      
      tlag <- u - 1
      if (tlag < 0)
        dV_lag <- 0
      else 
        dV_lag <- lagderiv(tlag, 3)  
      
      A[2] <- as.numeric(dV_lag > dV[2])
    }
    
    if( A[3] == 0 ) { 
      dV[3] <- G[3] - V[3]*m                              # vector net veg growth 
      
      tlag <- u - 1
      if (tlag < 0)
        dV_lag <- 0
      else 
        dV_lag <- lagderiv(tlag, 4)  
      
      A[3] <- as.numeric(dV_lag > dV[3])
    }
    
    dV <- G*(1 - A) - V*m                              # vector net veg growth 
    dX <- G*(A)
    
    dA <- A
    return( list( c(dR, dV, dA, dX ))) } )
}

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
r <- c(7, 4, 3.5) # max uptake rates mm of water per g of plant per day
K <- c(190, 35, 13)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.1                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.1         # rate of water evaporation and runoff mm per mm per day
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 

parms <- list( r = r, K = K, m =m, conversion = conversion, seedling_mass = seedling_mass, times = times)

R_init <- 200 
seeds_init <- c(2,2,2)

state <- c( R_init, seeds_init*parms$seedling_mass, 0, 0,0, 0, 0,0) 
times <- seq(0, 200, by = 0.01)

test <- dede( state, times = times, func = grow, parms = parms)
tail(test)

par(mfrow = c(1,1))
matplot( test[, 3:5], type = 'l')
matplot( test[, 3:5] + test[, 9:11], type = 'l')
matplot( test[, 9:11], type = 'l')

apply( test[, 9:11], 2, max)/(seeds_init*parms$seedling_mass)



