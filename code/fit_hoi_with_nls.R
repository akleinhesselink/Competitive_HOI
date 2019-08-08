rm(list = ls())

library(tidyverse) 

source('code/plotting_parameters.R')
source('code/simulation_functions.R')
source('code/phenomenological_models.R')

load( 'output/processed_results.rda')

fit4 <- fit3 <- fit2 <- fit1 <- list()

for( i in 1:3 ) { 
  
  temp_data <- sim_results %>% filter( n_comp < 2, species == !!paste0('Y', i))
  
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
    log(1/y) ~ log(model3(B1, B2, B3, parms = list(lambda = lambda, alpha = alpha, tau = tau))), 
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

# filter out two competitor runs where one 
# competitor is same as the focal species. 
temp_data <- 
  sim_results %>% 
  filter( n_comp < 3, (n_comp < 2) | ( get(paste0('B', str_extract(species, '\\d') )) < 1 )  )

fixed_pars3 <- lapply( fit3, get_fixed_pars)

hoi_init <- list( c(beta = 0.1, eta = 0.2), c( beta = 0.1, eta = 0.2), c( beta = -0.1, eta = 0.2))
hoi_lower <- list(c(0,0.1), c(0,0.1), c(-2,0.1))
hoi_upper <- list(c(2, 2), c(2, 2), c(-0.1, 2))

hoi_init2 <- lapply ( hoi_init, function(x) x[ 'beta' ] )
hoi_lower2 <- lapply ( hoi_lower, function(x) x[1] )
hoi_upper2 <- lapply ( hoi_upper, function(x) x[1 ] )



HOI_type1_fit <- HOI_type2_fit <- list()

for( i in 1:3) { 
  
  HOI_type1_fit[[i]] <- nls( 
    log( 1/y) ~  log( model_HOI_type_1(B1, B2, B3, h, beta = beta, eta = eta, fixed_parms = fixed_pars3[[i]] )) , 
    data = temp_data %>% filter( species == !!paste0( 'Y', i) & HOI == 1), 
    start = hoi_init[[i]], 
    lower = hoi_lower[[i]], 
    upper = hoi_upper[[i]], 
    algorithm = 'port', 
    control = list(maxiter = 1000, minFactor = 1/100000 )
  )
  
  HOI_type2_fit[[i]] <- 
    nls( 
      log( 1/y) ~  log( model_HOI_type_2[[i]](B1, B2, B3, beta = beta, fixed_parms = fixed_pars3[[i]] )) , 
      data = temp_data %>% filter( species == !!paste0( 'Y', i) & HOI == 1), 
      start = hoi_init2[[i]], 
      lower = hoi_lower2[[i]], 
      upper = hoi_upper2[[i]], 
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

pgrid$m1 <- NA
pgrid$m2 <- NA
pgrid$m3 <- NA
pgrid$m4 <- NA
pgrid$mHOI1 <- NA
pgrid$mHOI2 <- NA

for( i in 1:3 ) { 
  focal <- paste0('Y', i)
  pgrid$m1[pgrid$species == focal] <- 1/(exp( predict( fit1[[i]], newdata = pgrid[pgrid$species == focal, ]) ))
  pgrid$m2[pgrid$species == focal] <- 1/(exp( predict( fit2[[i]], newdata = pgrid[pgrid$species == focal, ]) ))
  pgrid$m3[pgrid$species == focal] <- 1/(exp( predict( fit3[[i]], newdata = pgrid[pgrid$species == focal, ]) ))
  pgrid$m4[pgrid$species == focal] <- 1/(exp( predict( fit4[[i]], newdata = pgrid[pgrid$species == focal, ]) ))
  
  pgrid$mHOI1[pgrid$species == focal] <- 1/(exp( predict( HOI_type1_fit[[i]], newdata = pgrid[pgrid$species == focal, ]) ))
  pgrid$mHOI2[pgrid$species == focal] <- 1/(exp( predict( HOI_type2_fit[[i]], newdata = pgrid[pgrid$species == focal, ]) ))
}



save( 
  temp_data, 
  pgrid,
  fit1, 
  fit2, 
  fit3, 
  fit4, 
  HOI_type1_fit, 
  HOI_type2_fit, file = 'output/model_fits.rda'
  )

