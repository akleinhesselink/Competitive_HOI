rm(list = ls())

library(tidyverse) 

source('code/plotting_parameters.R')
source('code/simulation_functions3.R')
source('code/phenomenological_models3.R')
load( 'output/processed_results3.rda')

nls.control(maxiter = 1000, tol = 1e-8, minFactor = (1/10)*(1/1024),
            printEval = FALSE, warnOnly = FALSE)

fit1pw <- fit2pw <- list()
fit2HOI <- fit1HOI <- list() 
fit2 <- fit1 <- list()

i <- 1
for( i in 1:3 ) { 
  # loop through focal species (i)
  
  temp_data <- 
    sim_results %>% 
    filter( n_comp < 2, species == !!paste0('Y', i))
  
  fit1pw[[i]] <- nls( 
    log(1/y) ~ log(model1(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits1, 
    lower = lowers1, 
    upper = uppers1, 
    algorithm = 'port'
  )

  fit2pw[[i]] <- nls( 
    log(1/y) ~ log(model2(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits2, 
    lower = lowers2, 
    upper = uppers2, 
    algorithm = 'port'
  )
}

for( i in 1:3) { 
  
  # fit HOI models 
  # loop through focal species (i)
  
  temp_data <- 
    sim_results %>% 
    filter( n_comp < 3, species == !!paste0('Y', i))
  
  fit1[[i]] <- nls( 
    log(1/y) ~ log(model1(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits1, 
    lower = lowers1, 
    upper = uppers1, 
    algorithm = 'port'
  )
  
  fit2[[i]] <- nls( 
    log(1/y) ~ log(model2(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits2, 
    lower = lowers2, 
    upper = uppers2, 
    algorithm = 'port'
  )
  
  fit1HOI[[i]] <- nls( 
    log(1/y) ~ log(model1_HOI(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, beta = beta, tau = tau))), 
    data = temp_data, 
    start = initsHOI, 
    lower = lowersHOI, 
    upper = uppersHOI,
    algorithm = 'port'
  )
  
}

fit1pw[[1]]
fit1pw[[2]]
fit1pw[[3]]

fit1[[1]]
fit1[[2]]
fit1[[3]]

fit2[[1]]
fit2[[2]]
fit2[[3]]

fit1HOI[[1]]
fit1HOI[[2]]
fit1HOI[[3]]

predicted <- pgrid

predicted$m1_HOI <- NA
predicted$m2 <- predicted$m1 <- NA
predicted$m1_pw <- NA

for( i in 1:3 ) { 
  focal <- paste0('Y', i)
  predicted$m1_pw[predicted$species == focal] <- 1/(exp( predict( fit1pw[[i]], newdata = predicted[predicted$species == focal, ]) ))
  #predicted$m2_pw[predicted$species == focal] <- 1/(exp( predict( fit2pw[[i]], newdata = predicted[predicted$species == focal, ]) ))
  
  predicted$m1[predicted$species == focal] <- 1/(exp( predict( fit1[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m2[predicted$species == focal] <- 1/(exp( predict( fit2[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m1_HOI[predicted$species == focal] <- 1/(exp( predict( fit1HOI[[i]], newdata = predicted[predicted$species == focal, ]) ))
  #predicted$m2_HOI[predicted$species == focal] <- 1/(exp( predict( fit2HOI[[i]], newdata = predicted[predicted$species == focal, ]) ))
}

save( 
  predicted,
  fit1pw, 
  fit2pw,
  fit1, 
  fit2, 
  fit1HOI,
  file = 'output/model_fits3.rda'
)

