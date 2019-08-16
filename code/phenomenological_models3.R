model1 <- function(B1, B2, B3, parms ){
  # Hassel Model 
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    
    (1/lambda)*(1 + t1 + t2 + t3 )^tau
  })
}

model2 <- function(B1, B2, B3, parms ){
  # Model 2, like model 1 but separate 
  # competitor terms are multiplied together
  # in the denominator 
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    
    (1/lambda)*((1 + t1)*(1 + t2)*(1 + t3))^tau
  })
}

model1_HOI <- function(B1, B2, B3, parms){ 
  # Hassel model with interspecific HOI terms 
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    HOI <- beta[1]*(B1*B2) + beta[2]*(B1*B3) + beta[3]*(B2*B3)
    (1/lambda)*(1 + t1 + t2 + t3 + HOI )^tau
  })
}

 
model2_HOI <- function(B1, B2, B3, parms){ 
  # Model 2 with HOI terms 
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    HOI1 <- beta[1]*(B1*B2)
    HOI2 <- beta[2]*(B1*B3)
    HOI3 <- beta[3]*(B2*B3)
    
    (1/lambda)*((1 + t1)*(1 + t2)*(1 + t3)*(1 + HOI1)*(1 + HOI2)*(1 + HOI3))^tau
  })
}


get_fixed_pars <- function(myfit ){ 
  
  x <- coef(myfit)
  list( 
    lambda = coef(myfit)[grep('lambda', names(coef(myfit)))], 
    alpha = coef(myfit)[grep('alpha', names(coef(myfit)))], 
    tau = coef(myfit)[ grep('tau', names(coef(myfit)))], 
    beta = coef(myfit)[ grep('beta', names(coef(myfit)))]
  )
}

# Initial values and parameter limits for nls fitting  
inits <- list( lambda = 10, alpha = rep(0.1, 3), tau = 0.9)
lowers <- unlist(inits)[] 
lowers[] <- c(1, rep(0, 3), 0) 
uppers <- lowers
uppers <- c(200, rep(10, 3), 5)

inits1 <- inits2 <- inits 
lowers1 <- lowers2 <- lowers
uppers1 <- uppers2 <- uppers

initsHOI <- inits1
initsHOI$beta <- rep(0, 3)
lowersHOI <- c( lowers1, beta = rep(-0.2, 3))
uppersHOI <- c( uppers, beta = rep(1, 3))