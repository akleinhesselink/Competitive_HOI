rm(list = ls())

library(deSolve)
library(tidyverse)

source('code/plotting_parameters.R')
source('code/phenomenological_models3.R')

load('output/model_fits4.rda')
load('output/simulate_coexistence.rda')

get_parms <- function(x) { 
  
  temp <- do.call(rbind, lapply( x, coef) )
  
  if( ncol(temp) == 5) {
    Lambda <- temp[, 'lambda']
    Alpha <- temp[, str_detect( colnames( temp  ) , 'alpha')]
    Tau <- temp[, 'tau']
    return( list(Lambda = Lambda, Alpha = Alpha, Tau = Tau)) 
  } else { 
    temp <- do.call(rbind, lapply( x, coef) )
    Lambda <- temp[, 'lambda']
    Alpha <- temp[, str_detect( colnames( temp  ) , 'alpha')]
    Beta <- temp[ , str_detect( colnames( temp  ) , 'beta' )]
    Tau <- temp[, 'tau']
    return( list(Lambda = Lambda, Alpha = Alpha, Beta = Beta , Tau = Tau))  
  }
}

reproduce <- function( N, parms){ 
  # N vector of abundances in year t-1 
  # parms are BH model parms
  N1 <- rep(NA, length(N))
  if(length(parms) == 3){ 
    with( parms , { 
      for( s in 1:length(N)){
        competitors <- N
        competitors[s] <- max( 0, N[s] - 1 ) # subtract one focal individual  
        N1[s] <- N[s]*(Lambda[s]/(1 + sum(Alpha[s,] * competitors ) )^(Tau[s]))
      }
      return( N1 )  
    })
  }else if( length( parms) == 4 ) { 
    with( parms , { 
      for( s in 1:length(N)){
        competitors <- N
        competitors[s] <- max( 0, N[s] - 1 ) # subtract one focal individual  
        HOI <- c( prod(competitors[1:2]), prod(competitors[c(1,3)]), prod(competitors[c(2,3)]), prod(competitors))
        if(ncol(Beta) == 3){ 
          N1[s] <- N[s]*(Lambda[s]/(1 + sum(Alpha[s,] * competitors ) + sum( Beta[s, ]*HOI[1:3] ) )^(Tau[s]))
        }else if(ncol(Beta) == 4){ 
          N1[s] <- N[s]*(Lambda[s]/(1 + sum(Alpha[s,] * competitors ) + sum( Beta[s, ]*HOI ) )^(Tau[s]))
        }
      }
      return( N1 )  
    })
  }
}


pw0_parms <- get_parms(fit1pw ) # fitted to single competitor gradients only 
pw1_parms <- get_parms( fit1) # fitted to multi-competitor gradients 
HOI1_parms <- get_parms(fit1HOI)
HOI2_parms <- get_parms(fit1HOI2)



simulated_comm <- lapply( community, function(x)  do.call( rbind, x )) 
phenom_comm <- simulated_comm
phenom_comm <- lapply( phenom_comm, function(x) { x[2:nrow(x), ] <- 0 ; return(x) } )

TotTime <- 200 
nSims <- length(phenom_comm)

phenom_comm_pw0 <- phenom_comm
phenom_comm_pw1 <- phenom_comm
phenom_comm_HOI1 <- phenom_comm
phenom_comm_HOI2 <- phenom_comm

for( j in 1:nSims ){ 
  for(t in 2:TotTime) { 
    phenom_comm_pw0[[j]][t, ] <- reproduce(phenom_comm_pw0[[j]][t-1, ], parms = pw0_parms)  
    phenom_comm_pw1[[j]][t, ] <- reproduce(phenom_comm_pw1[[j]][t-1, ], parms = pw1_parms)  
    phenom_comm_HOI1[[j]][t, ] <- reproduce(phenom_comm_pw1[[j]][t-1, ], parms = HOI1_parms)  
    phenom_comm_HOI2[[j]][t, ] <- reproduce(phenom_comm_pw1[[j]][t-1, ], parms = HOI2_parms)  
  }
}

phenom_comm_pw0
phenom_comm_pw1
phenom_comm_HOI1
phenom_comm_HOI2

compare_ts <- function(sim, phen) { 
  y_max <- max( sim, phen )
  
  matplot( sim , ylim = c(0, y_max), pch = 1)
  matplot( phen, type = 'l', add = T, lty = 1)
  legend( x = 150, y = y_max, legend = c('early', 'mid', 'late'), y.intersp = 0.2,
          pch = 1, 
          lty = 1, col = 1:3, cex = 1, bty = 'n')
}  

for( i in 1:length(simulated_comm)){ 
  compare_ts( simulated_comm[[i]], phenom_comm_pw0[[i]])
}
