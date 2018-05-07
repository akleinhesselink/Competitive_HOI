rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

parms_file <- 'data/mechanistic_parms.rda'

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
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)
save(parms, file = parms_file)
rm( list = names(parms) )

resource_curves <- plot_resource_uptake(parms)
ggsave(filename = 'figures/resource_uptake.png', resource_curves, height = 3, width = 4)

plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# plot time series and save 
seedlings <- c(1,0,1)
seedlings <- as.numeric(seedlings)
State <- c(parms$soil_m, seedlings*parms$seedling_mass)
out <- ode(y=State, times = seq( 1, parms$times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root)

ts_plot <- plot_timeseries(out)
ggsave( ts_plot, filename = 'figures/example_timeseries.png', height = 4, width = 4)

# -----------------------------------------------------

