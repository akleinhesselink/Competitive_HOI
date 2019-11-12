library(deSolve)

# Modeling Functions ------------------------------------------------------ #

g <- function(b, d, nu, p ) { 
  # Scales individual biomass "b" (in g) to 
  # root surface area (cm^2)  
  # d is root density (g/cm^3)  
  # nu scales root volume to root length, nu < 1 
  # p is proportion of total biomass allocated to roots

  SA <- ((p*b)/d)^nu
  SA[is.na(SA)] <- 0 
  SA
}



f <- function(R, Vmax, K){ 
  # Michaelis-Menton function for R uptake  
  # per unit SA 
  
  (Vmax*R)/(K + R) 
}


dBdu <- function(u, B, R, n, d, nu, p, Vmax, K, q, m) { 
  # B(u) is total biomass of population at time u
  # b(u) is the size per individual 
  
  b <- B/n
  n*(g(b, d, nu, p)*f(R, Vmax, K)*q - b*m)
} 

dRdu <- function(u, B, R, S, n, d, nu, p, Vmax, K) { 
  # R(u) is resource concentration at time u

  b <- B/n
  S - sum( n*g(b, d, nu, p)*f(R, Vmax, K)  )
} 


grow <- function(u, State, parms, ...){
  # Coupled differential equations for growth 
  # and resource uptake processes 
  with(parms , {
    R  <- State[1]                             
    B  <- State[2:length(State)]              

    dB <- dBdu(u, B, R, n, d, nu, p, Vmax, K, q, m)
    dR <- dRdu(u, B, R, S, n, d, nu, p, Vmax, K)

    return( list( c(dR, dB))) } )
}


root3 <- function(u, State, parms){ 
  # Resource uptake stops when resources are too 
  # low for positive net growth of individuals
  with(parms, 
       {
        b <- State[2:length(State)]/n
        
        out <- g(b, d, nu, p)*f(State[1], Vmax, K)*q - b*m
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




