###
#
# Fit models with three species HOI
#
##

rm(list = ls())

library(tidyverse) 

source('code/plotting_parameters.R')
source('code/simulation_functions3.R')
source('code/phenomenological_models3.R')
load( 'output/processed_results4.rda')

nls.control(maxiter = 1000, tol = 1e-8, minFactor = (1/10)*(1/1024),
            printEval = FALSE, warnOnly = FALSE)

fit1pw <- list()
fit1HOI2 <- fit1HOI <- list() 
fit1 <- list()

for( i in 1:3 ) { 
  # loop through focal species (i)
  temp_data <- 
    sim_results %>% 
    filter( n_comp < 2, species == !!paste0('Y', i))
  
  fit1pw[[i]] <- nls( 
    (1/y) ~ (model1(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits1, 
    lower = lowers1, 
    upper = uppers1, 
    algorithm = 'port'
  )
}

pw_coefs <- lapply( fit1pw, coef)
make_inits <- function(x){ 
  list( lambda = x['lambda'], alpha = x[2:4], tau = x['tau'])  
}

inits1 <- lapply(pw_coefs, make_inits)
initsHOI <- lapply( inits1, function(x) { x$beta <- c(0,0,0); return( x) } )
initsHOI2 <- lapply( inits1, function(x) { x$beta <- c(0,0,0,0); return( x) } )

i <- 1
for( i in 1:3) { 
  
  # fit HOI models 
  # loop through focal species (i)
  
  temp_data <- 
    sim_results %>% 
    filter(species == !!paste0('Y', i))
  uppers1

  if( i == 3) { 
    uppers1 <- c(10000, 50, 50, 50, 1) 
  }
  
  fit1[[i]] <- nls( 
    (1/y) ~ (model1(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits1[[i]], 
    lower = lowers1, 
    upper = uppers1, 
    algorithm = 'port'
  )
  
  fit1HOI[[i]] <- nls( 
    (1/y) ~ (model1_HOI(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, beta = beta, tau = tau))), 
    data = temp_data, 
    start = initsHOI[[i]], 
    lower = lowersHOI, 
    upper = uppersHOI,
    algorithm = 'port'
  )
  
  fit1HOI2[[i]] <- nls( 
    (1/y) ~ (model1_HOI2(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, beta = beta, tau = tau))), 
    data = temp_data, 
    start = initsHOI2[[i]], 
    lower = lowersHOI2, 
    upper = uppersHOI2,
    algorithm = 'port'
  )
  
}


fit1[[3]]
fit1HOI[[3]]
fit1HOI2[[3]]

fit1pw[[1]]
fit1pw[[2]]
fit1pw[[3]]

fit1[[1]]
fit1[[2]]
fit1[[3]]

fit1HOI[[1]]
fit1HOI[[2]]
fit1HOI[[3]]

fit1HOI2[[1]]
fit1HOI2[[2]]
fit1HOI2[[3]]

predicted <- pgrid

predicted$m1_HOI <- NA
predicted$m1_HOI2 <-NA
predicted$m1 <- NA
predicted$m1_pw <- NA

for( i in 1:3 ) { 
  focal <- paste0('Y', i)
  predicted$m1_pw[predicted$species == focal] <- 1/(exp( predict( fit1pw[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m1[predicted$species == focal] <- 1/(exp( predict( fit1[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m1_HOI[predicted$species == focal] <- 1/(exp( predict( fit1HOI[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m1_HOI2[predicted$species == focal] <- 1/(exp( predict( fit1HOI2[[i]], newdata = predicted[predicted$species == focal, ]) ))
}

save( 
  predicted,
  fit1pw, 
  fit1, 
  fit1HOI,
  fit1HOI2,
  file = 'output/model_fits4.rda'
)

