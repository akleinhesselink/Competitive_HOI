library(deSolve)

# Modeling Functions ------------------------------------------------------ #

K <- function(B, a = 1, k = 1 ){ 
  a/(2 - exp(-k*(B)))
}

f <- function(R, r, ...){ r*R/(K(...) + R) }                            # resource uptake rate. Saturates at r
dBdu <- function(u, B, R, r, q, m, a, k) { B*(q*f(R, r, B, a, k) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, a, k, n) { - sum(n*B*f(R, r, B, a, k) ) }  # resource 

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, r, q, m, a, k)
    dR <- dRdu(u, B, R, r, a, k, n)
    return( list( c(dR, dB))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K(State[2], a, k)/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K(State[2], a, k)/(q*r-m) < 0.000001) # logical vector of species to terminate
    State[2:length(State)][ terminate ] <- 0
    return(State)
  })
}
