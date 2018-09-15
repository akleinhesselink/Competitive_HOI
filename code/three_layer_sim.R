rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

# graphics themes ------------------------------------------------ # 

journal_theme <- 
  my_theme + 
  theme( axis.title = element_text(size = 12), 
         legend.text = element_text(size = 10), 
         legend.title = element_text(size = 12), 
         strip.text = element_text(size = 12), 
         axis.text = element_text(size = 10))

# Functions ------------------------------------------------------ #
f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r

dBdu <- function(u, B, R, r, w, K, q, m) { 
  
  B^(2/3)*q*sum(w*f(R, r, K)) - B*m
  
} # growth as a function of biomass and resource uptake

dRdu <- function(u, B, R, r, w, K) {
  
  - sum( B^(2/3)*w*f(R,r,K) )
  
} 

grow <- function(u, State, parms, ...){
  with(parms , {
    R1  <- State[1]                               # resource first
    R2  <- State[2]
    R3  <- State[3]
    B1  <- State[4]
    B2  <- State[5]
    B3  <- State[6]
    
    dB1 <- dBdu(u, B1, c(R1,R2,R3), r[1], w[1, ], K, q, m[1])
    dB2 <- dBdu(u, B2, c(R1,R2,R3), r[2], w[2, ], K, q, m[2])
    dB3 <- dBdu(u, B3, c(R1,R2,R3), r[3], w[3, ], K, q, m[3])
    
    dR1 <- dRdu(u, c(B1, B2, B3), R1, r, w[,1], K)
    dR2 <- dRdu(u, c(B1, B2, B3), R2, r, w[,2], K)
    dR3 <- dRdu(u, c(B1, B2, B3), R3, r, w[,3], K)
    
    return( list( c(dR1, dR2, dR3, dB1, dB2, dB3))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K/(q*r-m) < 0.000001) # logical vector of species to terminate
    State[2:length(State)][ terminate ] <- 0
    return(State)
  })
}

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
w <- matrix( c(1, 0, 0, 
               1, 1, 0, 
               1, 1, 1), 
             ncol = 3, 
             nrow = 3, 
             byrow = T)  # max uptake rates mm of water per g of plant per day

K <- 30
r <- c(3, 1, 0.8)
m <- c(0.05, 0.05, 0.05)                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
parms <- list( r = r, w = w, K = K, m = m, q = q, conversion = conversion, seedling_mass = seedling_mass, times = times)

R_init <- c(50, 200, 400) 
seeds_init <- c(1, 1, 1)

state <- c( R_init, seeds_init*seedling_mass)
test <- ode( state, times = 1:400, func = grow, parms = parms)

plot(test)
apply( test, 2, which.max)


plot(r, m)
plot(r, m, type= 'l')

out <- ode(state, times = 1:200, func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T), method = 'radau')

ts_plot <- plot_timeseries(out, sp_labs = c('Resource', '1', '2', '3'), mytheme = journal_theme + theme(legend.position = c(0.8, 0.5), axis.title = element_text(size = 24)))

ts_plot