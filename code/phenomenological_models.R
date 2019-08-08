# parameters for non-HOI models 
inits <- list( lambda = 10, alpha = rep(0.1, 3), tau = rep(0.1, 3))
lowers <- unlist(inits)[] 
lowers[] <- c(1, rep(0.01, 3), rep(0.001, 3)) 
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
    (1/lambda)*(1 + B1*alpha[1] + B2*alpha[2] + B3*alpha[3] )^tau
  })
}

model2 <- function(B1, B2, B3, parms ){
  # multiply separate species functions  
  with(parms, { 
    (1/lambda)*((1 + B1*alpha[1])*(1 + B2*alpha[2])*(1 + B3*alpha[3]))^tau
  })
}

model3 <- function(B1, B2, B3, parms ){
  with(parms, { 
    (1/lambda)*(1 + (B1*alpha[1])^tau[1] + (B2*alpha[2])^tau[2] + (B3*alpha[3])^tau[3] )
  })
}

model4 <- function(B1, B2, B3, parms ){
  with(parms, { 
    (1/lambda)*((1 + (B1*alpha[1])^tau[1])*(1 + (B2*alpha[2])^tau[2])*(1 + (B3*alpha[3])^tau[3]))^(-1)
  })
}

# HOI models 

model_HOI_type_1 <- function(B1, B2, B3, h, beta, eta, fixed_parms ){
  # formula is inverted 
  # fixed parms are from the pairwise fit 
  
  HOI <- beta*(h)^eta

  with(fixed_parms, {  
    (1/lambda)*(1 + (B1*alpha[1])^tau[1] + (B2*alpha[2])^tau[2] + (B3*alpha[3])^tau[3] + HOI)
  })
}

model_HOI_1 <- function(B1, B2, B3, beta, fixed_parms ){
  # formula is inverted 
  # fixed parms are from the pairwise fit 
  
  with(fixed_parms, {  
    HOI <- beta*( (B2*alpha[2])^tau[2]*(B3*alpha[3])^tau[3] )
    
    (1/lambda)*(1 + (B1*alpha[1])^tau[1] + (B2*alpha[2])^tau[2] + (B3*alpha[3])^tau[3] + HOI)
  })
}

model_HOI_2 <- function(B1, B2, B3, beta, fixed_parms ){
  # formula is inverted 
  # fixed parms are from the pairwise fit 
  
  with(fixed_parms, {  
    
    HOI <- beta*( (B1*alpha[1])^tau[1]*(B3*alpha[3])^tau[3] )
    (1/lambda)*(1 + (B1*alpha[1])^tau[1] + (B2*alpha[2])^tau[2] + (B3*alpha[3])^tau[3] + HOI)
    
  })
}



model_HOI_3 <- function(B1, B2, B3, beta, fixed_parms ){
  # formula is inverted 
  # fixed parms are from the pairwise fit 
  
  with(fixed_parms, {  
    
    HOI <- beta*( (B1*alpha[1])^tau[1]*(B2*alpha[2])^tau[2] )
    
    (1/lambda)*(1 + (B1*alpha[1])^tau[1] + (B2*alpha[2])^tau[2] + (B3*alpha[3])^tau[3] + HOI)
  })
}

model_HOI_type_2 <- list( model_HOI_1, model_HOI_2, model_HOI_3 )

get_fixed_pars <- function(myfit ){ 
  
  x <- coef(myfit)
  list( 
    lambda = coef(myfit)[grep('lambda', names(coef(myfit)))], 
    alpha = coef(myfit)[grep('alpha', names(coef(myfit)))], 
    tau = coef(myfit)[ grep('tau', names(coef(myfit)))] )
}




