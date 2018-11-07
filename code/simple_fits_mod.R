rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')

# Functions ------------------------------------------------------ #

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, r, K, q, m) { B*(q*f(R, r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K) { - sum(B*f(R,r, K)) } # resource (water)

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, r, K, q, m)
    dR <- dRdu(u, B, R, r, K)
    return( list( c(dR, dB))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K/(q*r-m) < 0.000001) # logical vector of species to terminate
    State[2:length(State)][ terminate ] <- 0
    return(State)
  })
}

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
r <- c(4.2, 2.6, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(150, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 200, length.out = 1000)
parms <- list( r = r, K = K, m =m, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)

save(parms,file = 'data/parms.rda')

resource_curves <- plot_resource_uptake(parms, spec_labs = species_labs, R = seq(0, 200, by = 0.01))

R <- 0:200
curves <- data.frame(R = R,  mapply(x = as.list(parms$r), y = as.list(parms$K), FUN = function(x, y) { f(R = R, x, y) }) )

curves <- 
  curves %>% 
  gather( species, uptake, starts_with('X')) %>% 
  mutate( species = factor(species, labels = species_labs))

resource_curves <- resource_curves + journal_theme + theme(legend.position = c(0.8, 0.3)) + ylab('Resource uptake rate')

ggsave(filename = 'figures//resource_uptake.png', resource_curves, height = 5.5, width = 6)

R_init <- 200 
seeds_init <- c(1,1,1)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass, seeds_init[3]*seedling_mass)

out <- ode(state, times = seq(1, 200, by = 0.01), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

ts_plots <- plot_timeseries(out, sp_labs = c('Resource', '1', '2', '3'), mytheme = journal_theme + theme(legend.position = c(0.8, 0.5), axis.title = element_text(size = 24)))

g1 <- ts_plots[[1]] + annotate(geom = 'text', 10, 220, label = 'a)', size = 5)
g2 <- ts_plots[[2]] + annotate(geom = 'text', 10, 2, label = 'b)', size = 5)

g3 <- resource_curves + 
  journal_theme  + 
  theme(legend.position = 'top', 
        legend.direction = 'vertical', 
        axis.text = element_blank(),
        axis.title = element_text(size = 24), 
        legend.title = element_text(size = 24), 
        legend.text = element_text(size = 20), 
        legend.key.size = unit(1.9, unit = 'lines'),
        legend.margin = margin(1.5,1,1,1, unit = 'lines')) + 
  ylab('Resource\nuptake rate') + 
  annotate( geom = 'text', 10, 2.4, label = 'c)', size = 5)

ts_plot <- 
  grid.arrange( 
  g1, g2, g3,
  widths = c(1,1), 
  heights = c(1,1.2),
  layout_matrix = rbind( c(1,3), 
                         c(2,3))
  )

ggsave( ts_plot, filename = 'figures//example_timeseries.png', height = 5, width = 7)

# Run response surface experiments --------------------------- # 

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
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  state[4] <- B_init[i,3]*parms$seedling_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, 200, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  
}

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

save(sim_results, file = 'data/sim_results.rda')
