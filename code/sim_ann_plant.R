library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)
rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

load('data/experiments.rda')
load('data/out.rda')
load('data/parms.rda')
load('data/seedling_mass.rda')
load('data/conversion.rda')

load('data/result_fit.rda')
load('data/experiments.rda')
load('data/mech_eq.rda')

ann_plant_mod <- function(x, form, pars) { 
  lambda <- pars[1]   
  alphas <- pars[2:(length(pars)-1)]
  tau <- pars[length(pars)]
  mm <- model.matrix(as.formula( form ), x )
  y <- lambda*(1 + mm%*%alphas)^(tau)    
  return(y)
}

results[[1]][[2]]
results[[2]][[2]]
results[[3]][[2]]

State <- c(1,0,1)
#as.matrix( apply( experiments, 2, median))
State <- data.frame( t(State) ) 
names( State) <- c('N1', 'N2', 'N3')

form1 <- paste('~ -1 + N1 + N2 + N3')
form2 <- paste(form1, '+ I(N1*N2) + I(N1*N3) + I(N2*N3)')
all_forms <- c(form1, form2)

nspp <- length(State)
time <- 100

out1 <- out2 <- data.frame( matrix(NA, ncol = nspp, nrow = time))
out1[1, ] <- out2[1, ] <- as.numeric( State)
names(out1) <- names(out2) <- c('N1', 'N2', 'N3')

pars1 <- lapply( results, function(x) x[[1]]$par )
pars2 <- lapply( results, function(x) x[[2]]$par )

for(j in 2:time) { 
  for( i in 1:nspp){ 
    out1[j, i] <- out1[j-1, i]*ann_plant_mod(out1[j-1, ], form = form1, pars = pars1[[i]])
    out2[j, i] <- out2[j-1, i]*ann_plant_mod(out2[j-1, ], form = form2, pars = pars2[[i]])
  }
} 

ap1_eq <- out1[nrow(out1), ]
ap2_eq <- out2[nrow(out2), ]
mech_eq

ap2_eq
ap1_eq

(ap1_eq - mech_eq)/mech_eq
(ap2_eq - mech_eq)/mech_eq

#
State <- as.numeric( ap1_eq)
State[2] <- 1
#as.matrix( apply( experiments, 2, median))
State <- data.frame( t(State) ) 
names( State) <- c('N1', 'N2', 'N3')

form1 <- paste('~ -1 + N1 + N2 + N3')
form2 <- paste(form1, '+ I(N1*N2) + I(N1*N3) + I(N2*N3)')
all_forms <- c(form1, form2)

nspp <- length(State)
time <- 1000

out1 <- out2 <- data.frame( matrix(NA, ncol = nspp, nrow = time))
out1[1, ] <- out2[1, ] <- as.numeric( State)
names(out1) <- names(out2) <- c('N1', 'N2', 'N3')

pars1 <- lapply( results, function(x) x[[1]]$par )
pars2 <- lapply( results, function(x) x[[2]]$par )

for(j in 2:time) { 
  for( i in 1:nspp){ 
    out1[j, i] <- out1[j-1, i]*ann_plant_mod(out1[j-1, ], form = form1, pars = pars1[[i]])
    out2[j, i] <- out2[j-1, i]*ann_plant_mod(out2[j-1, ], form = form2, pars = pars2[[i]])
  }
} 

ap1_eq
mech_eq
out1[time, ]
matplot(out1)

#
State <- as.numeric( ap1_eq)
State[3] <- 1
#as.matrix( apply( experiments, 2, median))
State <- data.frame( t(State) ) 
names( State) <- c('N1', 'N2', 'N3')

form1 <- paste('~ -1 + N1 + N2 + N3')
form2 <- paste(form1, '+ I(N1*N2) + I(N1*N3) + I(N2*N3)')
all_forms <- c(form1, form2)

nspp <- length(State)
time <- 1000

out1 <- out2 <- data.frame( matrix(NA, ncol = nspp, nrow = time))
out1[1, ] <- out2[1, ] <- as.numeric( State)
names(out1) <- names(out2) <- c('N1', 'N2', 'N3')

pars1 <- lapply( results, function(x) x[[1]]$par )
pars2 <- lapply( results, function(x) x[[2]]$par )

for(j in 2:time) { 
  for( i in 1:nspp){ 
    out1[j, i] <- out1[j-1, i]*ann_plant_mod(out1[j-1, ], form = form1, pars = pars1[[i]])
    out2[j, i] <- out2[j-1, i]*ann_plant_mod(out2[j-1, ], form = form2, pars = pars2[[i]])
  }
} 

ap1_eq
mech_eq
out1[time, ]
matplot(out1)


