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
times <- 200             # length of simulation in days 
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
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)
rm( list = names(parms) )

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# -------- simulate multiple years -------------------------- # 
t <- 500 # number of years 

experiments <- expand.grid( c(0, 10), c(0, 10), c(0, 10))
experiments <- split(experiments, 1:nrow(experiments))
experiments <- lapply( experiments, as.numeric)
results <- lapply( experiments, run_multi_gen, t = t, parms = parms, tol = 1e-5)

eqs <- do.call( rbind, lapply( results, function(x) { cc <- x[complete.cases(x), ]; cc[nrow(cc), ] } ) )
experiments <- do.call(rbind, experiments)

results <- as.data.frame (do.call(rbind, results))

population_sims <- 
  tibble::rownames_to_column(results, var = 'id') %>% 
  separate(id, c('experiment', 'generation'), sep = '\\.')



save(parms, file = 'data/parms.rda')
save(eqs, file = 'data/eqs.rda')
save(population_sims, file = 'data/population_sims.rda')
save(experiments, file = 'data/experiments.rda')




