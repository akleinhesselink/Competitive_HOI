library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(parallel)

rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 125             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.9, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(98, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedlings <- c(1, 0, 0)      # number of seedlings 
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q)
rm( list = names(parms) )

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

State <- c(soil_m, seedlings*seedling_mass)
out <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )

plot_timeseries(out, parms, col = my_colors)

seeds <- c(1,1,1)
State <- c(soil_m, seeds*seedling_mass)
out <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
plot_timeseries(out, parms, col = my_colors)

# -------- simulate multiple years -------------------------- # 
t = 200 # number of years 
seedlings <- c(0.5,0.5,0.5)
State <- c(soil_m, seedlings*seedling_mass)
use <- out <- list(NA)
fecundity <- data.frame(N1 = rep(NA, t), N2 = rep(NA, t), N3 = rep(NA, t))

for( i in 2:t){ 
  out[[i-1]] <- ode(y=State, times = seq(1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
  use[[i-1]] <- matrix(NA, nrow = nrow(out[[i-1]]), ncol = 3)
  for( k in 1:3){ 
    use[[i-1]][, k] <- out[[i-1]][,2 + k]*f(out[[i-1]][, 2], parms$r[k], parms$K[k])
  }
  max_biomass  <- apply( out[[i-1]][, c(3:5)], 2, max)
  fecundity[i-1,] <- (max_biomass*conversion)/seedling_mass
  State <- as.numeric(c(soil_m, fecundity[i-1, ]*seedling_mass))
}

par(mfrow = c(1,1))
matplot(fecundity, type = 'l', col = my_colors)

# fit parameters: 




