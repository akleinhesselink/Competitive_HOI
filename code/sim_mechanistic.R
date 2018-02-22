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

# -------- simulate annual plant experiments -------------------------- # 
experiments <- expand.grid(c(0, seq(1, 25, 1)), c(0, seq(1, 35, 1)), c(0, seq(1, 40, 1))) # response surface experiment 
experiments <- experiments[-1, ]
experiments <- experiments %>% distinct()

experiments <- 
  rbind( experiments [ rowSums( experiments[, 1:3] == 0) == 2 , ] , 
       experiments [ apply( experiments, 1, function(x) any( x == 1 ) ) & rowSums( experiments[, 1:3] == 0) == 1 , ], 
       experiments [ experiments[, 1] > 15 & experiments[, 2] > 25 & experiments[, 3] > 30, ])


#experiments <- experiments[ rowSums(experiments) < 30, ] 

run_experiment <- function(seedlings, soil_m, seedling_mass) { 
  seedlings <- as.numeric(seedlings)
  State <- c(soil_m, seedlings*seedling_mass)
  ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
}

start1 <- proc.time()
out <- mclapply( 
  split(experiments, f = 1:nrow(experiments)), 
  run_experiment, 
  soil_m, 
  seedling_mass, mc.cores = 4
  )
t1 <- proc.time() - start1
t1

save(out, file = 'data/out.rda')
save(experiments, file = 'data/experiments.rda')
save(parms, file = 'data/parms.rda')
save(seedling_mass, file = 'data/seedling_mass.rda')
save(conversion, file = 'data/conversion.rda')
