rm(list = ls())
nspp <- 2
times <- 100 
A <- matrix(c(0, 5.1, 
              0.2, 0), 2, 2, byrow = T)

N_init <- c(1, 0)
N <- matrix(NA, times, nspp, byrow = T)
N[1, ] <- N_init

for( i in 2:times){ 
  Nint <- A %*% N[i-1, ]
  N[i, ] <- A %*% Nint 
}

matplot( N, type = 'l' )

# 

g <- function(x, b = 1){ 
  exp(-b*x)
}

logistic <- function(x){ 
   1/(1 + exp(-x)) 
}

survival <- function(x, beta0 = 2, beta1 = -0.5){ 
  
  logistic(beta0 + beta1*x) 
  
}


curve(logistic, 0, 3)
beta0 <- 2
beta1 <- -3

curve( survival(x), 0, 10) 
survival 
x*(max_fecundity - alpha_f %*% x)



