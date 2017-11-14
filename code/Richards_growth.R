lgrowth <- function(t, B, parms){ 
  with (parms, { 
    out = r*B*(1 - alpha%*%B)    
    return(list(out))
  })
}



library(deSolve)

nspp <- 3 
r <- c(0.15, 0.12, 0.11)
K <- 1/r*2

alpha <- matrix(mean(K), nrow = nspp, ncol = nspp)
diag(alpha) <- K
alpha <- 1/alpha

B <- rep(0.01, nspp)

parms <- list(alpha = alpha, r = r)

out <- ode(B, times = 1:300, func = lgrowth, parms = parms )

plot(out)

