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
tiny <- .Machine$double.eps
times <- 125             # length of simulation in days 
soil_m <- 500            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
alpha <- 0.1            # factors determining trade-off curve 
beta <- 6
gamma <- 0
delta <- 6
cc <- 10
TO_pars <- list( alpha = alpha, beta = beta , gamma  = gamma, delta = delta, cc = cc)

r <- c(7, 4) # max uptake rates mm of water per g of plant per day
m <- 0.02                # tissue respiration and loss rate g per g per day 
q <- 0.1                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.1         # rate of water evaporation and runoff mm per mm per day
seedlings <- c(1, 1)      # number of seedlings 
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 100, length.out = 1000)

parms <- list( r = r, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q)
parms$K <- TO_fun(parms$r, TO_pars)
plot_transpiration(parms, my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms,  my_colors)

State <- c(soil_m, seedlings*seedling_mass)
out <- ode(y=State, func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root, times = seq( 1, times, 0.1))
plot_timeseries(out, parms, col = my_colors)

# simulate_ESS 
par(mfrow = c(1,1))

my_seq <- seq(0.4, 10, length.out = 30)

r <- my_seq

start <- proc.time()
N_hat <- mclapply( r, 
        function(x) { 
          parms$r <- x;
          find_N_hat(fun = seed_production, 
                     init_size = seedling_mass, 
                     R_state = soil_m, 
                     parms = parms,
                     TO_fun = TO_fun, 
                     TO_pars = TO_pars, 
                     func = grow,  
                     events = list(func = event, root = TRUE), 
                     rootfun = root, 
                     times = seq( 1, times, 0.1), method = 'adams')}, mc.cores = 4)
proc.time() - start 

test_parms <- parms
test_parms$r <- 3

fitness(1e-3, fun  = seed_production, init_size = seedling_mass, 
         R_state = soil_m, 
         parms = test_parms,
         TO_fun = TO_fun, 
         TO_pars = TO_pars, 
         func = grow,  
         events = list(func = event, root = TRUE), 
         rootfun = root, 
         times = seq( 1, times, 0.1), method = 'adams')

find_N_hat_old(fun = seed_production, 
           init_size = seedling_mass, 
           R_state = soil_m, 
           parms = test_parms,
           TO_fun = TO_fun, 
           TO_pars = TO_pars, 
           func = grow,  
           events = list(func = event, root = TRUE), 
           rootfun = root, 
           times = seq(1, times, 0.1), method = 'adams',
           my_range = c(1e-7, 100))

proc.time() - start 

start = proc.time()
find_N_hat(fun = seed_production, 
            init_size = seedling_mass, 
            R_state = soil_m, 
            parms = test_parms,
            TO_fun = TO_fun, 
            TO_pars = TO_pars, 
            func = grow,  
            events = list(func = event, root = TRUE), 
            rootfun = root, 
            times = seq( 1, times, 0.1), method = 'adams',
            my_range = c(1e-7, 100))
proc.time() - start 



N_hat <- unlist(lapply(N_hat, function(x) x$minimum))

sims <- expand.grid(r1 = my_seq, r2 = my_seq)

sims$N_hat <- N_hat 

out <- mclapply(split(sims, 1:nrow(sims)), FUN = function(x, ... ){parms$r <- unlist(x[1:2]); fitness(x = unlist(c(x[3], tiny)), parms = parms, fun = seed_production, ... )}, 
         init_size = seedling_mass, 
         R_state = soil_m, 
         func = grow,
         TO_fun = TO_fun, 
         TO_pars = TO_pars, 
         events = list(func = event, root = TRUE), 
         rootfun = root, 
         times = seq(1, times, 0.1), mc.cores = 4)

out <- do.call(rbind, out)
sims$lambda2 <- out[,2]

ggplot(sims, aes(x = r1, y = r2, z = lambda2, fill = lambda2 > 0)) + 
  geom_tile() + geom_contour() + geom_abline(aes(slope = 1, intercept = 0))

#
ggplot( sims %>% filter( r1 < 7.5 & r1 > 7.3 ), aes( x = r2, y = lambda2 ) ) + geom_line() + geom_hline(aes( yintercept = 0), linetype = 2)


plot( sims[ sims$r1 < 7.5 & sims$r1 > 7.3, ]$lambda2 , type = 'l')


unique(sims$r1)

fitness(c(0,tiny),
        parms = parms, 
        fun = seed_production, 
        init_size = seedling_mass, 
        R_state = soil_m, 
        func = grow,
        TO_fun = TO_fun, 
        TO_pars = TO_pars, 
        events = list(func = event, root = TRUE), 
        rootfun = root, 
        times = seq(1, times, 0.1))
