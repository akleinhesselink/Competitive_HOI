#
# Run simulations for Appendix A 
# 

rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
library(gridExtra)
library(grid)

source('code/simulation_functions.R')
source('code/plotting_parameters.R')
source('code/phenomenological_models.R')
source('code/process_sim_data.R')

load('output/parms.rda')

nsims <- 5 
nsps  <- 3
max_d <- 36

# new funct ions to adjust the parameters controlling resource uptake 
curve_parms <- list( a = 0.05, b = 2000, r = 0.6)

trade_off <- function(m, curve_parms ) { 
  with(curve_parms, { 
    r*(dgamma(m, shape = a, scale = b))
  })
}

m3s <- seq( parms$m[3], parms$m[2] - 0.005, length.out = nsims)
m1s <- seq( parms$m[1], parms$m[2] + 0.005, length.out = nsims)
m2s <- parms$m[2]

m_parms <- cbind( m1= m1s, m3 = m3s)

trade_off_trials <- 
  data.frame( trial = 1:nsims, m1s, m2s, m3s) %>% 
  gather( species, m, starts_with('m')) %>%
  mutate( species = factor( species, labels = species_labs) ) %>%
  mutate( d = trade_off(m, curve_parms )) 

save(curve_parms, trade_off, trade_off_trials, file = 'output/trade_off_trials.rda')

# ---------- set up experimental gradient ------------------- # 
max_d <- max_d # set above 
min_d <- 0 
dgrad <- sort( c(2,3,(seq(sqrt(min_d), sqrt(max_d), 1))^2) )

B_init <- expand.grid( 
  B1 = dgrad, 
  B2 = dgrad,
  B3 = dgrad)

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 ) %>%  # filter out three species cases 
  filter( B1 + B2 + B3 < 50) # filter out cases where total density is high  

state <- c( R = parms$R0, B1 = 0, B2 = 0, B3 = 0)
res <- list()

for( k in 1:nsims){
  
  parms$m <- as.numeric( c( m_parms[k, 1], parms$m[2], m_parms[k, 2]))
  parms$d <- as.numeric( trade_off(parms$m, curve_parms))
  out <- list()
  for( i in 1:nrow(B_init)){
    
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

  sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
  sim_results <- data.frame(B_init, sim_results )

  res[[k]] <- 
    sim_results %>% 
    rename( "X1" = R, "X2" = B1.1, "X3" = B2.1, "X4" = B3.1) %>% 
    mutate( trial = k) 
}

# Process the data from each trial -------- # 

processed_data <- lapply( res, process_results, parms = parms)

# Fit models to each scenario ------------- # 

fits <- list()

initsHOI <- list(lambda = 10, alpha = rep(0.1, 3), tau = 0.5, beta = rep(0, 3))  
lowersHOI <- c(lambda = 1, alpha1 = -1e-2, alpha2 = -1e-2, alpha3 = -1e-3, tau = 1e-04, beta = -1e-02, beta2 = -1e-02, beta3 = -1e-02)

for( k in 1:nsims){ 
  
  sim_results <- processed_data[[k]]
  fitHOI <- list()
  for( i in 1:nsps) { 
    
    # fit HOI models 
    # loop through focal species (i)
    
    temp_data <- 
      sim_results %>% 
      filter( n_comp < 3, species == !!paste0('Y', i))
    
    fitHOI[[i]] <- nls( 
      (1/y) ~ (model1_HOI(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, beta = beta, tau = tau))), 
      data = temp_data, 
      start = initsHOI, 
      lower = lowersHOI, 
      upper = uppersHOI,
      algorithm = 'port'
    )
  }
  
  fits[[k]] <- fitHOI
}

res <- do.call( bind_rows, 
         lapply( fits, lapply, function(x) data.frame( t(coef(x)))) )

res$species <- factor(1:3,  labels = species_labs, ordered = T)
res$trial <- sort(rep( c(1:nsims), nsps))

tradeoff_results <- 
  res %>% 
  gather( par, value, lambda:beta3) %>% 
  mutate( par_type = str_extract(par, '[a-z]+')) %>% 
  group_by( trial, species, par_type ) %>% 
  mutate( par_avg = mean(abs(value)))


save(tradeoff_results, file = 'output/trade_off_results.rda')


