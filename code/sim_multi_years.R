library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)

rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 125             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.3, 2.9, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(110, 35, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
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
t = 100 # number of years 
comp_grad <- c(0:1)  # number of competitors
experiments <- expand.grid(as.list(rep(list ( comp_grad), 3))) # response surface experiment 
experiments <- experiments[-1, ]
experiments <- experiments[order(rowSums(experiments)), ]

fecundity <- phenology <- data.frame( matrix( NA, nrow = t, ncol = 3))
per_capita_use <- use <- out <- list()
pb <- txtProgressBar(min = 1, nrow(experiments)*t -1, style = 3)

counter <- 1
for(j in 1:nrow(experiments)){ 
  seedlings <- as.numeric(experiments[j, ])
  for( i in 2:t){ 
    setTxtProgressBar(pb , counter)
    State <- c(soil_m, seedlings*seedling_mass)
    out[[i-1]] <- ode(y=State, times = seq(1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
    use[[i-1]] <- matrix(NA, nrow = nrow(out[[i-1]]), ncol = 3)
    for( j in 1:3){ 
      use[[i-1]][, j] <- out[[i-1]][,2 + j]*f(out[[i-1]][, 2], parms$r[j], parms$k[j])
    }
    phenology[i-1,] <- apply( out[[i-1]][, c(3:5)], 2, find_phenology)
    max_biomass  <- apply( out[[i-1]][, c(3:5)], 2, max)
    seedlings <- fecundity[i-1,] <- (max_biomass*conversion)/seedling_mass
    counter = 1 + counter
  }
  matplot(fecundity, type = 'l')
}




