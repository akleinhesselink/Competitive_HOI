library(deSolve)

# Modeling Functions ------------------------------------------------------ #

g <- function(b, b0, d, nu){    
  # Scales biomass to surface area (SA)
  # of resource acqusition organs (roots & leaves). 
  # SA scales sub-linearly with biomass. 
  # Species differ in tissue density (d), 
  # species with greater lower density tissues 
  # have greater SA. 
  
  SA <- b0*(b/d)^(nu) 
  SA[is.na(SA)] <- 0 
  SA 
}

f <- function(R, Vmax, K){ 
  # Michaelis-Menton function for R uptake  
  # per unit SA 
  
  (Vmax*R)/(K + R) 
}


dBdu <- function(u, B, R, n, b0, d, nu, Vmax, K, q, m) { 
  # B(u) is total biomass of population at time u
  # b(u) is the size per individual 
  
  b <- B/n
  n*(g(b, b0, d, nu)*f(R, Vmax, K)*q - b*m)
} 

dRdu <- function(u, B, R, n, b0, d, nu, Vmax, K) { 
  # R(u) is resource concentration at time u

  b <- B/n
  - sum( n*g(B/n, b0, d, nu)*f(R, Vmax, K)  )
} 


grow <- function(u, State, parms, ...){
  # Coupled differential equations for growth 
  # and resource uptake processes 
  with(parms , {
    R  <- State[1]                             
    B  <- State[2:length(State)]              

    dB <- dBdu(u, B, R, n, b0, d, nu, Vmax, K, q, m)
    dR <- dRdu(u, B, R, n, b0, d, nu, Vmax, K)

    return( list( c(dR, dB))) } )
}


root3 <- function(u, State, parms){ 
  # Resource uptake stops when resources are too 
  # low for positive net growth of individuals
  with(parms, 
       {
        b <- State[2:length(State)]/n
        
        out <- g(b, b0, d, nu)*f(State[1], Vmax, K)*q - b*m
        out[ out == 0 ]  <- 1
        
        return(out)
       }
  )
}


event <- function(u, State, parms) {
  # Terminate species (set B(u) = 0) once 
  # resources are too low for positive net growth
  
  terminate <- root3(1, State, parms) < 1e-9  # logical vector of species to terminate
  State[2:length(State)][ terminate ] <- 0
  return(State)
}




