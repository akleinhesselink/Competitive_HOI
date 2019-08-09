rm(list = ls())

library(tidyverse) 

source('code/plotting_parameters.R')
source('code/simulation_functions.R')
source('code/phenomenological_models.R')
load( 'output/processed_results.rda')

fit4 <- fit3 <- fit2 <- fit1 <- list()

for( i in 1:3 ) { 
  
  temp_data <- 
    sim_results %>% 
    filter( n_comp < 2, species == !!paste0('Y', i))
  
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
  fit3[[i]] <- nls(
    log(1/y) ~ log(model3(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits3, 
    lower = lowers3, 
    upper = uppers3, 
    algorithm = 'port'
  )
  fit4[[i]] <- nls(
    log(1/y) ~ log(model4(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
    data = temp_data, 
    start = inits4, 
    lower = lowers4, 
    upper = uppers4, 
    algorithm = 'port'
  )
  
}

fit1[[1]]
fit1[[2]]
fit1[[3]]

fit2[[1]]
fit2[[2]]
fit2[[3]]

fit3[[1]] 
fit3[[2]]
fit3[[3]]

fit4[[1]]
fit4[[2]]
fit4[[3]]

fixed_pars3 <- lapply( fit3, get_fixed_pars)


HOI_type1_fit <- HOI_type2_fit <- list()
for( i in 1:3) { 
  
  temp_data <- 
    sim_results %>% 
    filter( n_comp < 3, species == !!paste0('Y', i))
  
  HOI_type1_fit[[i]] <- nls( 
    log( 1/y) ~  log( model_HOI_type_1(B1, B2, B3, beta = beta, eta = eta, fixed_parms = fixed_pars3[[i]] )) , 
    data = temp_data,
    start = hoi_init[[i]], 
    lower = unlist( hoi_lower[[i]]), 
    upper = unlist( hoi_upper[[i]]), 
    algorithm = 'port'
  )
  
  HOI_type2_fit[[i]] <- 
    nls( 
      log( 1/y) ~  log( model_HOI_type_2(B1, B2, B3, beta = beta, fixed_parms = fixed_pars3[[i]] )) , 
      data = temp_data,
      start = hoi_init[[i]][1], 
      lower = hoi_lower[[i]][1:3], 
      upper = hoi_upper[[i]][1:3], 
      algorithm = 'port', 
      control = list(maxiter = 1000, minFactor = 1/100000 )
    )
  
}

HOI_type1_fit[[1]]
HOI_type1_fit[[2]]
HOI_type1_fit[[3]]

HOI_type2_fit[[1]]
HOI_type2_fit[[2]]
HOI_type2_fit[[3]]

predicted <- pgrid

predicted$m1 <- NA
predicted$m2 <- NA
predicted$m3 <- NA
predicted$m4 <- NA
predicted$mHOI1 <- NA
predicted$mHOI2 <- NA

for( i in 1:3 ) { 
  focal <- paste0('Y', i)
  predicted$m1[predicted$species == focal] <- 1/(exp( predict( fit1[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m2[predicted$species == focal] <- 1/(exp( predict( fit2[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m3[predicted$species == focal] <- 1/(exp( predict( fit3[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$m4[predicted$species == focal] <- 1/(exp( predict( fit4[[i]], newdata = predicted[predicted$species == focal, ]) ))
  
  predicted$mHOI1[predicted$species == focal] <- 1/(exp( predict( HOI_type1_fit[[i]], newdata = predicted[predicted$species == focal, ]) ))
  predicted$mHOI2[predicted$species == focal] <- 1/(exp( predict( HOI_type2_fit[[i]], newdata = predicted[predicted$species == focal, ]) ))
}

save( 
  predicted,
  fit1, 
  fit2, 
  fit3, 
  fit4, 
  HOI_type1_fit, 
  HOI_type2_fit, file = 'output/model_fits.rda'
  )

