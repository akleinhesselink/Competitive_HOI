rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, S, r, K, q, m) { B*(q*f(R + S[u], r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K) { - sum(B*f(R,r, K)) } # resource (water)

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, S, r, K, q, m)
    dR <- dRdu(u, B, R, r, K)
    return( list( c(dR, dB))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K/(q*r-m) < 0) # logical vector of species to terminate
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
R <- seq(0, 500, length.out = 1000)
S <- rep(0, times)

parms <- list( r = r, K = K, m =m , q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times, S = S)
R_init <- 200 

experiments <- 
  expand.grid(B1 = c(0,1), B2 = c(0,1), B3 = c(0,1), pulse = 1:times) %>% 
  filter( B1 + B2 + B3 == 1) %>% 
  arrange( B1, B2, B3, pulse )

out <- list()

for(i in 1:nrow(experiments)){ 
  delta <- 0.01
  seeds_init <- as.numeric(experiments[i, 1:3])
  S <- rep(0, parms$times)
  S[experiments$pulse[i]] <- delta
  parms$S <- S

  state <- c(R_init, 
             seeds_init[1]*seedling_mass, 
             seeds_init[2]*seedling_mass, 
             seeds_init[3]*seedling_mass)

  out[[i]] <- ode(state, 
             times = 1:200,
             func = grow, 
             parms = parms, 
             rootfun = root, 
             event = list(func = event, root = T), 
            method = 'radau')

}

baseline_experiments <- 
  expand.grid(B1 = c(0,1), B2 = c(0,1), B3 = c(0,1)) %>% 
  filter( B1 + B2 + B3 == 1) %>% 
  arrange( B1, B2, B3)

baseline_out <- list()
parms$S <- rep(0, parms$times)

for(i in 1:nrow(baseline_experiments)){ 
  
  seeds_init <- as.numeric(baseline_experiments[i, 1:3])
  
  state <- c(R_init, 
             seeds_init[1]*seedling_mass, 
             seeds_init[2]*seedling_mass, 
             seeds_init[3]*seedling_mass)
  
  baseline_out[[i]] <- ode(state, 
                  times = 1:200,
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), 
                  method = 'radau')
  
}

final_biomass <- do.call(rbind, lapply( out, function(x) apply( x , 2, max)[3:5] ))
baseline_biomass <- do.call(rbind, lapply(baseline_out, function(x) apply( x, 2, max)[3:5]))

sens1 <- (final_biomass[ rowSums(final_biomass[, c(2,3) ]) == 0 , 1] - baseline_biomass[3, 1])/delta
sens2 <- (final_biomass[ rowSums(final_biomass[, c(1,3) ]) == 0 , 2] - baseline_biomass[2, 2])/delta
sens3 <- (final_biomass[ rowSums(final_biomass[, c(1,2) ]) == 0 , 3] - baseline_biomass[1, 3])/delta

par(mfrow = c(1,1))
plot(sens2, type = 'l', col = my_colors[2])
points(sens3, type = 'l', col = my_colors[3])
points(sens1, type = 'l')

baseline_out[[3]]
imp3 <- f(R = baseline_out[[1]][,2], r = parms$r[3], K = parms$K[3] )*(baseline_out[[1]][, 5] > 0 )
imp2 <- f(R = baseline_out[[2]][,2], r = parms$r[2], K = parms$K[2] )*(baseline_out[[2]][, 4] > 0 )
imp1 <- f(R = baseline_out[[3]][,2], r = parms$r[1], K = parms$K[1] )*(baseline_out[[3]][, 3] > 0 )

plot(imp1, type = 'l')
points(imp2, type = 'l', col = my_colors[2])
points(imp3, type = 'l', col = my_colors[3])
# 

plot( sens1*imp3, type = 'l', col = my_colors[3])
points( sens1*imp2, type = 'l', col = my_colors[2])
points( sens1*imp1, type = 'l', col = my_colors[1])

plot( sens2*imp2, type = 'l', col = my_colors[2])
points( sens2*imp1, type = 'l', col = my_colors[1])
points( sens2*imp3, type = 'l', col = my_colors[3])

plot( sens3*imp3, type = 'l', col = my_colors[3])
points( sens3*imp2, type = 'l', col = my_colors[2])
points( sens3*imp1, type = 'l', col = my_colors[1])

plot(sens2*imp1, type = 'l', col = my_colors[2])
points(sens3*imp1, type = 'l', col = my_colors[1])

barplot( cbind( mean(sens2[imp1 > 0]), mean(sens3[imp1 > 0]) ) )
barplot( cbind( mean(sens3[imp1 > 0]), mean(sens3[imp2 > 0])))

barplot( cbind( mean(sens2[imp1 > 0]), mean(sens2[imp3 > 0])))
barplot( cbind( mean(sens3[imp1 > 0]), mean(sens2[imp1 > 0])))
barplot( cbind( mean(sens3[imp1 > 0]), mean(sens2[imp1 > 0])))


plot( ((imp1 > 0)*sens3)/max(((imp1 > 0)*sens3)), type = 'l')
points( imp1/max(imp1), type = 'l')

plot( ((imp1 > 0)*sens2)/max(((imp1 > 0)*sens2)), type = 'l')
points( imp1/max(imp1), type = 'l')


plot(sens2*imp3, type = 'l', col = my_colors[2])
points(sens1*imp3, type = 'l', col = my_colors[1])

plot(sens1*imp2, type = 'l', col = my_colors[2])
points(sens1*imp3, type = 'l', col = my_colors[3])


# Compare species 2 and 3 interactions with and without species one 

species_one_experiments <- 
  expand.grid(B1 = c(0,1), B2 = 1, B3 = 1, pulse = 1:parms$times) %>% 
  arrange( B1, B2, B3)

species_one_experiments

species_one_out <- list()
parms$S <- rep(0, parms$times)
R_init <- 200 

for(i in 1:nrow(species_one_experiments)){ 
  
  seeds_init <- as.numeric(species_one_experiments[i, 1:3])
  
  S <- rep(0, parms$times)
  S[experiments$pulse[i]] <- 1
  parms$S <- S
  
  state <- c(R_init, 
             seeds_init[1]*seedling_mass, 
             seeds_init[2]*seedling_mass, 
             seeds_init[3]*seedling_mass)
  
  species_one_out[[i]] <- ode(state, 
                           times = 1:200, 
                           func = grow, 
                           parms = parms, 
                           rootfun = root, 
                           event = list(func = event, root = T), 
                           method = 'radau')
  
}

species_one_baseline <- 
  expand.grid(B1 = c(0,1), B2 = 1, B3 = 1) %>% 
  arrange( B1, B2, B3)

species_one_baseline

species_one_baseline_out <- list()
parms$S <- rep(0, parms$times)

for(i in 1:nrow(species_one_baseline)){ 
  
  seeds_init <- as.numeric(species_one_baseline[i, 1:3])
  
  state <- c(R_init, 
             seeds_init[1]*seedling_mass, 
             seeds_init[2]*seedling_mass, 
             seeds_init[3]*seedling_mass)
  
  species_one_baseline_out[[i]] <- ode(state, 
                           times = 1:200, 
                           func = grow, 
                           parms = parms, 
                           rootfun = root, 
                           event = list(func = event, root = T), 
                           method = 'radau')
  
}


species_one_biomass <- do.call(rbind, lapply( species_one_out, function(x) apply( x , 2, max)[3:5] ))
species_one_biomass[1:200, 3]

species_one_baseline <- do.call(rbind, lapply( species_one_baseline_out, function(x) apply( x , 2, max)[3:5] ))

par(mfrow = c(1,1))
plot( species_one_biomass[1:200, 3] , type = 'l')

sens3_0 <- species_one_biomass[ 1:200, 3] - species_one_baseline[1, 3]
plot( sens3_0 , type = 'l')

sens2_0 <- species_one_biomass[ 1:200, 2] - species_one_baseline[1, 2]
plot( sens2_0 , type = 'l')

sens3_1 <- species_one_biomass[201:400, 3] - species_one_baseline[2, 3]
plot( sens3_1 , type = 'l')

sens2_1 <- species_one_biomass[201:400, 2] - species_one_baseline[2, 2]
plot( sens2_1 , type = 'l')


par(mfrow = c(1,1))
plot(sens3_0, type = 'l', col = my_colors[3], lty = 2)
points(sens2_0, type  = 'l', col = my_colors[2], lty = 2)

points(sens2_1, type = 'l', col = my_colors[2])
points(sens3_1, type = 'l', col = my_colors[3])


imp_3_0 <- (species_one_baseline_out[[1]][, 5])*f( R = species_one_baseline_out[[1]][,2], r = parms$r[3], K = parms$K[3] )
plot(imp_3_0, type = 'l')

imp_2_0 <- (species_one_baseline_out[[1]][, 4])*f( R = species_one_baseline_out[[1]][,2], r = parms$r[2], K = parms$K[2] )
plot(imp_2_0, type = 'l')

imp_3_1 <- (species_one_baseline_out[[2]][, 5])*f( R = species_one_baseline_out[[2]][,2], r = parms$r[3], K = parms$K[3] )
plot(imp_3_1, type = 'l')

imp_2_1 <- (species_one_baseline_out[[2]][, 4])*f( R = species_one_baseline_out[[2]][,2], r = parms$r[2], K = parms$K[2] )
plot(imp_2_1, type = 'l')




plot(imp_2_0, type = 'l', col = my_colors[2])
points(imp_3_0, type = 'l', col = my_colors[3])
points(imp_2_1, type = 'l', col = my_colors[2], lty = 2)
points(imp_3_1, type = 'l', col = my_colors[3], lty = 2)

plot(sens3_0*imp_2_0, type = 'l', ylim = c(-2e-04, 6e-04))
points(sens3_1*imp_2_1, type = 'l', lty = 2)

plot(sens2_0*imp_3_0, type = 'l')
points(sens2_1*imp_3_1, type = 'l', lty = 2)
