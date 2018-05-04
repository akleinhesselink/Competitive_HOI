rm(list = ls())

dBdt <- function(B, R, q, u, m) { 
  B*(q*u*R - m )
}

dRdt <- function(R, B1, B2, u){ 
  
  0 - B1*u[1]*R - B2*u[2]*R  
  
}

grow <- function(t, state, parms ){ 
  
  R <- state[1]
  B1 <- state[2]
  B2 <- state[3]
  
  with(parms, { 
    dB1 <- dBdt( B1, R, q, u[1], m[1])
    dB2 <- dBdt( B2, R, q, u[2], m[2])
    dR  <- dRdt( R, B1, B2, u)
    
    
    return(list(c(dR = dR, dB1 = dB1, dB2 = dB2)))
  })
}
  
root <- function(t, state, parms) with(parms, { c(state[1]*q*u[1] - m[1], state[1]*q*u[2] - m[2]) } )

event <- function(t, state, parms) {
  with(parms, {
    B1_term <- state[1]*q*u[1] - m[1] < 1e-7 
    B2_term <- state[1]*q*u[2] - m[2] < 1e-7
    
    
    if(B1_term){ 
      state[2] <- 0 
    }
    
    if(B2_term){ 
      state[3] <- 0 
    }
  
      return(state)
  })
}

library(deSolve)

R_init <- 500
B_init <- c(1, 2)

parms <- list(u = c(0.00101, 0.001), 
               m = c(0.089, 0.001), q = 0.5)

state <- c(R_init, B_init)

out <- ode(state, times = 1:200, func = grow, parms = parms )

plot(out)

out <- ode(state, times = 1:300, func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

plot( out )
out

state <- c(R_init, c(3.5,0))
out <- ode(state, times = seq(1,200, length = 500), func = grow, parms = parms, method = 'radau',
           rootfun = root, event = list(func = event, root = T))

plot( out )

B_init <- expand.grid( B1 = c(0, seq(1, 20, by = 0.5)), B2 = c(0, seq(1, 20, by = 0.5)))
B_init <- B_init[-1,]
B_init
out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]
  state[3] <- B_init[i,2]
  
  out[[i]] <- ode(state, times = seq(1,200, length = 500), func = grow, parms = parms, 
             rootfun = root, event = list(func = event, root = T), method = 'radau')
}

plot(out[[1]])
plot(out[[21]])

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))

results <- data.frame(B_init, results )

plot( results$X2[1:19 ], type = 'l' )

library(tidyverse)
library(stringr)

results %>% 
  filter( B1 == 0 ) %>% 
  ggplot( aes( x = B2, y = X3)) + 
  geom_point() + 
  geom_abline(aes(intercept = 243.5, slope = 1 ))


results <- 
  results %>% 
  mutate( Y1 = X2/B1, Y2 = X3/B2) %>% 
  select( - X1)


results <- 
  results %>% 
  gather( species, y, Y1, Y2)  %>% 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) 


results <- results %>% 
  filter( !is.na(y))

results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))


form1 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2)^tau.'
form2 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.*I(B1*B2))^tau.'
form3 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.[1]*I(B1*B2) + beta.[2]*I(B2^2))^tau.'
formB <- 'y ~ lambda/(1 + B1^alpha.[1] + B2^alpha.[2])^tau.'

init_pars1 <- list( alpha. = c(1, 1), tau. = 1)
init_pars2 <- list( alpha. = c(1, 1), beta. = 0, tau. = 1)
init_pars3 <- list( alpha. = c(1, 1), beta. = c(0,0), tau. = 1)
init_parsB <- list( alpha. = c(1, 1), tau. = 1)

upper1 <- c(10, 10, 2)
lower1 <- c(0, 0, 0)
upper_HOI_1 <- c(10, 10, 3, 2)
lower_HOI_1 <- c(0, 0, 0, 0, 0)
upper_HOI_2 <- c(10, 10, 3, 3, 2)
lower_HOI_2 <- c(0, 0, 0, 0, 0, 0)

fit_1 <- nls(form1, data = results %>% filter( species == 'Y1'), start = init_pars1, algorithm = 'port', lower = lower1, upper = upper1)
fit_2 <- nls(form1, data = results %>% filter( species == 'Y2'), start = init_pars1, algorithm = 'port', lower = lower1, upper = upper1)

fit_1
fit_2

fit_1_HOI <- nls(form2, data = results %>% filter( species == 'Y1'), start = init_pars2, algorithm = 'port', lower = lower_HOI_1, upper = upper_HOI_1)
fit_2_HOI <- nls(form2, data = results %>% filter( species == 'Y2'), start = init_pars2, algorithm = 'port', lower = lower_HOI_1, upper = upper_HOI_1)

fit_1_HOI
fit_2_HOI

fit_1_HOI2 <- nls(form3, 
                  data = results %>% filter( species == 'Y1'), 
                  start = init_pars3, 
                  algorithm = 'port', 
                  lower = lower_HOI_2, 
                  upper = upper_HOI_2)

fit_2_HOI2 <- nls(form3, 
                  data = results %>% filter( species == 'Y2'), 
                  start = init_pars3, 
                  algorithm = 'port', 
                  lower = lower_HOI_2, 
                  upper = upper_HOI_2)

fit_1_B <- nls(formB, data = results %>% filter( species == 'Y1'), start = init_parsB, algorithm = 'port', lower = c(0, 0, 0), upper = c(5, 5, 2))
fit_2_B <- nls(formB, data = results %>% filter( species == 'Y2'), start = init_parsB, algorithm = 'port', lower = c(0, 0, 0), upper = c(5, 5, 2))

fit_2
fit_2_HOI
fit_2_HOI2

fit_1
fit_1_HOI
fit_1_HOI2

results <- 
  results %>% 
  mutate( pred1 = ifelse( species == 'Y1', predict( fit_1), predict( fit_2))) %>% 
  mutate( pred_HOI = ifelse( species == 'Y1', predict( fit_1_HOI), predict( fit_2_HOI2)))

p1 <- results %>% 
  group_by( species ) %>% 
  filter( B2 %in% c(0, 1, 2, 3, 4)) %>% 
  do(gg = ggplot( data = ., aes(x = B1, y = y, color = factor(B2) )) + 
       geom_point() + 
       geom_line(aes(y = pred1), linetype = 2) + 
       geom_line(aes(y = pred_HOI), linetype = 3) + 
       ylab( unique(.$species)))


p1$gg[[1]]
p1$gg[[2]]

fit_2
fit_2_HOI


p2 <- results %>% 
  group_by( species ) %>% 
  do(gg = ggplot( data = ., aes(x = B1, y = B2, z = lambda/y )) + 
       stat_contour() )


p2$gg[[1]]
p2$gg[[2]]


cond <- expand.grid(B = seq(0, 400, by = 50), R = seq(0, 400, by = 50))
parms

cond <- cond %>% 
  mutate(dBdt = dBdt(B, R, 0.5, 0.001, 0.08), 
         dRdt = dRdt(R, B1 = B, B2 = 0, u = c(0.001, 0)))


cond %>% 
  ggplot( aes( x = B, y = R)) + 
  geom_segment(aes(x = B, xend = B + dBdt, y = R, yend = R + dRdt), arrow = arrow(length = unit(3, 'pt')))

