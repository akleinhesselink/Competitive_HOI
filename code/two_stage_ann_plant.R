rm(list = ls())
logit <- function( t ) { 
    
  1/ ( 1 + exp(-t))
  
}

inv.logit <- function( p ) { 
  
  log( p / ( 1 - p )  )
}

logit(0.5)
inv.logit( 0.5 )

a11 <- 0.1
a12 <- 0.05
a21 <- 0.05
a22 <- 0.1 

n1 <- 10
n2 <- 20

survival <- function( n, alpha ) { 
  logit(1 - alpha[1]*n[1] - alpha[2]*n[2] - alpha[3]*n[3] )
}


fecundity <- function(a, lambda, beta ) { 
  
  lambda*(1 + a[1]*beta[1] + a[2]*beta[2] + a[3]*beta[3])^-1 
}

alpha <- matrix( c(0.01, 0.1, 0, 
                   0.0, 0.01, 0, 
                   0.0, 0.0, 0.01), 3,3,byrow = T )

beta <- matrix( c(0.01, 0, 0.0, 
                   0.0, 0.01, 1, 
                   0.0, 0.0, 0.01), 3,3,byrow = T )

alpha
beta 

one_gen <- function( a, alpha, beta, lambda ){   
  n <- NA
  n[1] <- a[1]*fecundity( a, lambda[1], beta[1, ] )
  n[2] <- a[2]*fecundity( a, lambda[2], beta[2, ] )
  n[3] <- a[3]*fecundity( a, lambda[3], beta[3, ] )
  
  a[1] <- n[1]*survival( n, alpha[1,])
  a[2] <- n[2]*survival( n, alpha[2,])
  a[3] <- n[3]*survival( n, alpha[3,])  
  
  return(a)
}

rs <- expand.grid( 0:10, 0:10, 0:10 )
out <- rs
out[] <- NA

for( i in 1:nrow(rs)) { 
  a <- as.numeric(rs[i, ])
  out[i, ] <- one_gen(a = a, alpha, beta, lambda = c(10, 10, 10))
}

y <- out/rs 

y <- y + rnorm( length(y), 0, 0.1)

testform <- 'y ~ lambda.*(1 + alpha.[1]*x1 + alpha.[2]*x2 + alpha.[3]*x3)^-1'  
testform2 <- 'y ~ lambda.*(1 + alpha.[1]*x1 + alpha.[2]*x2 + alpha.[3]*x3 + beta.[1]*x2*x3)^-1'

mydat <- data.frame( rs[ , 1:3])

names(mydat) <- c('x1', 'x2', 'x3')

sp <- 1
mydat$y <- y[, sp]
mydat <- mydat[ mydat[, sp] > 0 , ] 
mydat[, sp] <- mydat[,sp] - 1 

testfit <- nls(testform, data = mydat, start = list(lambda. = 10, alpha. = c(0,0,0)))
testfit

testfit <- nls( testform2, data = mydat, start = list(lambda. = 10, alpha. = c(0,0,0), beta. = c(0.2)))
testfit


