# parameters for non-HOI models 
inits <- list( lambda = 10, alpha = rep(0.1, 3), tau = rep(0.1, 3))
lowers <- unlist(inits)[] 
lowers[] <- c(1, rep(0, 3), rep(0.001, 3)) 
uppers <- lowers
uppers <- c(200, rep(10, 3), rep(5, 3))

inits1 <- inits2 <- inits3 <- inits4 <- inits 
lowers1 <- lowers2 <- lowers3 <- lowers4 <- lowers
uppers1 <- uppers2 <- uppers3 <- uppers4 <- uppers

inits2$tau <- inits1$tau <- inits1$tau[1]
lowers2 <- lowers1 <- c( lowers1[1:4], tau = 0.001)
uppers2 <- uppers1 <- c( uppers1[1:4], tau = 5)

# simple pairwise models 
model1 <- function(B1, B2, B3, parms ){
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    
    (1/lambda)*(1 + t1 + t2 + t3 )^tau
  })
}

model2 <- function(B1, B2, B3, parms ){
  # multiply separate species functions  
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    
    (1/lambda)*((1 + t1)*(1 + t2)*(1 + t3))^tau
  })
}

model3 <- function(B1, B2, B3, parms ){
  with(parms, { 
    t1 <- (B1*alpha[1])^tau[1] 
    t2 <- (B2*alpha[2])^tau[2]
    t3 <- (B3*alpha[3])^tau[3]
    
    (1/lambda)*(1 + t1 + t2 + t3 )
  })
}

model4 <- function(B1, B2, B3, parms ){
  with(parms, { 
    t1 <- (B1*alpha[1])^tau[1] 
    t2 <- (B2*alpha[2])^tau[2]
    t3 <- (B3*alpha[3])^tau[3]
    
    (1/lambda)*((1 + t1)*(1 + t2)*(1 + t3))
  })
}

model5 <- function(B1, B2, B3, parms ){ 
  
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    
    lambda*(1 - t1 - t2 - t3)
  })
}

# HOI models 
model_HOI_type_00 <- function(B1, B2, B3, parms){ 
  with(parms, { 
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    HOI <- beta[1]*(B1*B2) + beta[2]*(B1*B3) + beta[3]*(B2*B3)
    (1/lambda)*(1 + t1 + t2 + t3 + HOI )^tau
  })
}


model_HOI_type_0 <- function(B1, B2, B3, beta, fixed_parms ) { 
  with(fixed_parms, {  
    t1 <- (B1*alpha[1]) 
    t2 <- (B2*alpha[2])
    t3 <- (B3*alpha[3])
    
    HOI <- beta[1]*(B1*B2) + beta[2]*(B1*B3) + beta[3]*(B2*B3)
    
    (1/lambda)*(1 + t1 + t2 + t3 + HOI)^tau
  })
}

model_HOI_type_1 <- function(B1, B2, B3, beta, eta, fixed_parms ){
  # formula is inverted 
  # fixed parms are from the pairwise fit 
  
  with(fixed_parms, {  
    t1 <- (B1*alpha[1])^tau[1] 
    t2 <- (B2*alpha[2])^tau[2]
    t3 <- (B3*alpha[3])^tau[3]
    
    HOI <- beta[1]*(B1*B2)^eta[1] + beta[2]*(B1*B3)^eta[2] + beta[3]*(B2*B3)^eta[3]
    
    (1/lambda)*(1 + t1 + t2 + t3 + HOI)
  })
}

model_HOI_type_2 <- function(B1, B2, B3, beta, fixed_parms ){
  # formula is inverted 
  # fixed parms are from the pairwise fit 
  
  with(fixed_parms, {  
    t1 <- (B1*alpha[1])^tau[1] 
    t2 <- (B2*alpha[2])^tau[2]
    t3 <- (B3*alpha[3])^tau[3]
    
    HOI <- beta[1]*(t1*t2) + beta[2]*(t1*t3) + beta[3]*(t2*t3)
    
    (1/lambda)*(1 + t1 + t2 + t3 + HOI)
  })
}

get_fixed_pars <- function(myfit ){ 
  
  x <- coef(myfit)
  list( 
    lambda = coef(myfit)[grep('lambda', names(coef(myfit)))], 
    alpha = coef(myfit)[grep('alpha', names(coef(myfit)))], 
    tau = coef(myfit)[ grep('tau', names(coef(myfit)))] )
}


beta_init <- matrix( rep(0.1, 9), 3, 3)
beta_init <- beta_init*matrix( c(rep(1, 6), -1, 1,1), 3,3, byrow = T )

eta_init <- matrix( rep(0.5, 9), 3,3)
beta_lower <- (beta_init - abs(beta_init))*10
beta_upper <- (beta_init + abs(beta_init))*10

eta_lower <- eta_init*0.1
eta_upper <- eta_init*2

hoi_init <- list(list(beta = beta_init[1, ], eta = eta_init[1, ]),
                 list(beta = beta_init[2, ], eta = eta_init[2, ]),
                 list(beta = beta_init[3, ], eta = eta_init[3, ]))

hoi_upper <- list( c(beta_upper[1,], eta_upper[1, ]), 
                   c(beta_upper[2,], eta_upper[2, ]), 
                   c(beta_upper[3,], eta_upper[3, ]))

hoi_lower <- list( c(beta_lower[1,], eta_lower[1, ]), 
                   c(beta_lower[2,], eta_lower[2, ]), 
                   c(beta_lower[3,], eta_lower[3, ]))


