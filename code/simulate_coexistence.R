rm(list = ls())

library(deSolve)
library(tidyverse)
source('code/simulation_functions3.R')
source('code/process_sim_data3.R')

load( 'output/parms3.rda')

get_equilibrium <- function( init_density = c(B1 = 1, B2 = 0, B3 = 0), TotTime = 100, ...  ){ 

  out <- list()
  out[[1]] <- init_density 
  
  for( t in 2:TotTime){
    parms$n <- as.numeric( out[[t-1]][1:3] )
    parms$n[parms$n == 0 ] <- 1
  
    state <- c(R = parms$R0, out[[t-1]]*parms$seed_mass)
  
    res <- ode(state, 
                  times = seq(1, parms$U, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root3, 
                  event = list(func = event, root = T))
  

    out[[t]] <- apply( res, 2, max)[3:5]*parms$conversion/parms$seed_mass
  }
  
  return( out )
}

init_df <- expand.grid( B1 = c(0,1), B2 = c(0,1), B3 = c(0,1)) %>% 
  filter( (B1 + B2 + B3) == 1 )

mono <- list()

for( j in 1:nrow(init_df)){ 
  mono[[j]] <- get_equilibrium(
      init_density = unlist(init_df[j,]), 
      TotTime = 200, 
      parms = parms
      )
} 

end_density <- unlist( lapply( mono, tail, 1 ) , recursive = F)
end_density <- do.call( rbind, end_density )

equilibs <- end_density[ end_density != 0  ] 

init_invade <- 
  expand.grid( B1 = c(0,1,equilibs[1]), B2 = c(0, 1, equilibs[2]), B3 = c(0,1, equilibs[3]) ) %>% 
  filter( (B1 > 1) + (B2 > 1) + (B3 > 1 ) == 1 ) %>% 
  filter( (B1 < 1) + (B2 < 1) + (B3 < 1 ) == 1 ) 

invasions <- list()

for( j in 1:nrow(init_invade)){ 
  invasions[[j]] <- get_equilibrium(
    init_density = unlist(init_invade[j,]), 
    TotTime = 50, 
    parms = parms
  )
} 


invasion_result <- do.call( rbind, lapply( lapply( invasions, function(x) do.call(rbind, x )) , tail, 1 ))

# simulate invasion into all possible 2-species community 

init_invade2 <- invasion_result
init_invade2[init_invade2 == 0 ] <- 1

invasions2 <- list()

for( j in 1:nrow(init_invade2)){ 
  invasions2[[j]] <- get_equilibrium(
    init_density = unlist(init_invade2[j,]), 
    TotTime = 50, 
    parms = parms
  )
} 

invasion_results2 <- do.call(rbind , lapply(  lapply( invasions2, function(x) do.call(rbind, x)) , tail, 1 ) )

save( mono, invasion_result, invasion_results2, file = 'output/simulate_coexistence.rda')
