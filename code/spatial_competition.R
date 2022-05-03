rm(list = ls())
nls.control(maxiter = 1000, tol = 1e-8, minFactor = (1/10)*(1/1024),
            printEval = FALSE, warnOnly = FALSE)

# Parameters
p <- 0.5  # fraction habitat 1 
K <- 2 # different kinds of sites 
S <- 3 # species diversity 
t <- 100
area <- 100

# habitat 1 interactions
A1 <- matrix( c(1, 0.9,   0.1, 
               0,   1, 0.9, 
               0.9,   0, 1), S, S, byrow = T)

lambda1 <- c(5,1,4)

# habitat 2 interactions 
A2 <- matrix( c(1, 0.9,   0.1, 
                0,   1, 0.9, 
                0.9,   0, 1), S, S, byrow = T)
lambda2 <- c(1,5,4)

BH <- function( x , pars ){ 
  with(pars, {
    (x*(lambda/( 1 + A %*% x)))*area
  })
}

pars1 <- list( lambda = lambda1, A = A1, area = area)
pars2 <- list( lambda = lambda2, A = A2, area = area)

N2 <- N1 <- matrix( rep(NA, t*S) , t, S)
N1[1, ] <- c(0.01, 0.01, 0.01)
N2 <- N1 
N <- N1 + N2 



for( i in 1:(t-1)) { 
  
  temp1 <- rpois( 3, BH(N1[i, ], pars1) )/area
  temp2 <- rpois( 3, BH(N2[i, ], pars2) )/area
  N[i + 1, ] <- temp1 + temp2 
  N1[i + 1, ] <- p*N[i + 1, ]
  N2[i + 1, ] <- (1-p)*N[i + 1, ]
}

matplot(1:t,  N[1:t, ], type = 'l')
legend('topright', legend = c(1,2,3), lty = c(1,2,3), col = c(1,2,3))


x <- N[-t, ]
y <- N[-1, ]/x

pop_dat <- data.frame( x = N[-t,], y = x/(N[-1, ]))

model1 <- function(x.1, x.2, x.3, parms ){
  # Hassel Model 
  with(parms, { 
    t1 <- (x.1*alpha[1]) 
    t2 <- (x.2*alpha[2])
    t3 <- (x.3*alpha[3])
    
    (1/lambda)*(1 + t1 + t2 + t3 )
  })
}

model2 <- function(x.1, x.2, x.3, parms ){
  # Hassel Model 
  with(parms, { 
    t1 <- (x.1*alpha[1]) 
    t2 <- (x.2*alpha[2])
    t3 <- (x.3*alpha[3])
    
    h12 <- (x.1*x.2*beta[1])
    h13 <- (x.1*x.3*beta[2])
    h23 <- (x.2*x.3*beta[3])

    (1/lambda)*(1 + t1 + t2 + t3 + h12 + h13 + h23 )
  })
}

lowers <- as.matrix( expand.grid(0, c(0, -1), c(0,-1), c(0,-1)))
m1 <- list() 
for( i in 1:nrow(lowers)){ 
  
  m1[[i]] <- 
    try( 
      nls( 
        log(y.1) ~ log(model1( x.1, x.2, x.3, parms = list( lambda = lambda, alpha = alpha))), 
        data = pop_dat, 
        start = list( lambda = 1, alpha = c(1, 1, 1)), 
        lower = lowers[i, ], 
        upper = c(10, 5, 5, 5), 
        algorithm = 'port' 
        ), 
      silent = T)
}  

m1 <- m1[ ! lapply(m1, class) == 'try-error']
m1.top <- m1[ which.min( unlist( lapply( m1, function(x) x$m$deviance()) )) ] 
lowers <- expand.grid( 0, c(0, -1), c(0,-1),c(0,-1), c(0,-0.5), c(0,-0.5), c(0,-0.5))
lowers
m2 <- list() 
for( i in 1:nrow(lowers)){ 
  
  m2[[i]] <- 
    try( 
      nls( 
        log(y.1) ~ log(model2( x.1, x.2, x.3, parms = list( lambda = lambda, alpha = alpha, beta = beta))), 
        data = pop_dat, 
        start = list( lambda = 1, alpha = c(1, 1, 1), beta = c(0,0,0)), 
        lower = lowers[i, ], 
        upper = c(10, 5, 5, 5, 1, 1, 1), 
        algorithm = 'port' 
      ), 
      silent = T)
}  


m2 <- m2[ ! lapply(m2, class) == 'try-error']
m2.top <- m2[ which.min( unlist( lapply( m2, function(x) x$m$deviance()) )) ] 

m1.top
m2.top
