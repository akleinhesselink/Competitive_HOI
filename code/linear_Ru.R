rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)

source('code/figure_pars.R')

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
    dB1 <- dBdt( B1, R, q[1], u[1], m[1])
    dB2 <- dBdt( B2, R, q[2], u[2], m[2])
    dR  <- dRdt( R, B1, B2, u)
    
    
    return(list(c(dR = dR, dB1 = dB1, dB2 = dB2)))
  })
}
  
root <- function(t, state, parms) with(parms, { c(state[1]*q[1]*u[1] - m[1], state[1]*q[2]*u[2] - m[2]) } )

event <- function(t, state, parms) {
  with(parms, {
    B1_term <- state[1]*q[1]*u[1] - m[1] < 1e-7 
    B2_term <- state[1]*q[2]*u[2] - m[2] < 1e-7
    
    
    if(B1_term){ 
      state[2] <- 0 
    }
    
    if(B2_term){ 
      state[3] <- 0 
    }
  
      return(state)
  })
}

par(mfrow = c(1,1))
R_init <- 450
B_init <- c(1, 1)

parms <- list(u = c(0.001, 0.001), 
              m = c(0.06, 0.001), 
              q = c(0.6, 0.42))


state <- c(R_init, B_init)

out <- ode(state, times = 1:100, func = grow, parms = parms )
matplot(out[, -c(1:2)], type = 'l', lty = c(1,1))

out <- ode(state, times = seq(1, 100, by = 0.1), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T), method = 'radau')
matplot(out[, -c(1:2)], type = 'l', lty = c(1,1))

df <- data.frame(out)  
names(df) <- c('time', 'R', '1', '2')

df <- 
  df %>% 
  gather( species, biomass, `1`:`2`) 

ts_plot <- 
  df %>% 
  ggplot( aes( x = time, y = biomass, color = species )) + 
  geom_line() + 
  scale_color_manual(values = my_colors) + 
  my_theme + 
  theme(axis.text = element_blank())

ts_plot

B_init <- expand.grid( B1 = c(0, seq(1, 10, by = 1)), B2 = c(0, seq(1, 10, by = 1)))
B_init <- expand.grid( B1 = c(0, seq(1, 20000, by = 1000)), B2 = 0 )
B_init <- B_init[-1,]

out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]
  state[3] <- B_init[i,2]
  
  out[[i]] <- ode(state, times = seq(1, 100, by = 0.01), func = grow, parms = parms, 
             rootfun = root, event = list(func = event, root = T), method = 'radau')
}
B_init

test_case <- apply( out[[1]], 2, max)
test_case
plot(out[[1]])

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))


plot(B_init$B1, results[,3])
abline(0,1)

results <- data.frame(B_init, results )




results <- 
  results %>% 
  mutate( Y1 = X2/B1, Y2 = X3/B2) %>% 
  select( - X1)



results <- 
  results %>% 
  gather( species, y, Y1, Y2)  %>% 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0))

results <- results %>% 
  filter( !is.na(y))

results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

form1 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2)^tau.'
form2 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.*I(B2^2))^tau.'
form3 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.*I(B1*B2))^tau.'
form4 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.[1]*I(B2^2) + beta.[2]*I(B1*B2))^tau.'

init_pars1 <- list( alpha. = c(1, 1), tau. = 1)
init_pars2 <- list( alpha. = c(1, 1), beta. = 0, tau. = 1)
init_pars3 <- list( alpha. = c(1, 1), beta. = 0, tau. = 1)
init_pars4 <- list( alpha. = c(1, 1), beta. = c(0,0), tau. = 1)

upper1 <- c(10, 10, 2)
lower1 <- c(0, 0, 0)
upper_HOI_1 <- c(10, 10, 10, 2)
lower_HOI_1 <- c(0, 0, 0, 0, 0)
upper_HOI_2 <- c(10, 10, 10, 10, 2)
lower_HOI_2 <- c(0, 0, 0, 0, 0, 0)


fit_1 <- nls(form1, 
             data = results %>% filter( species == 'Y1'), 
             start = init_pars1, 
             algorithm = 'port', 
             lower = lower1, 
             upper = upper1)

fit_2 <- nls(form1, 
             data = results %>% filter( species == 'Y2'), 
             start = init_pars1, 
             algorithm = 'port', 
             lower = lower1, 
             upper = upper1)


fit_1_HOI_1 <- nls(form2, 
                 data = results %>% filter( species == 'Y1'), 
                 start = init_pars2, 
                 algorithm = 'port', 
                 lower = lower_HOI_1, 
                 upper = upper_HOI_1)

fit_2_HOI_1 <- nls(form2, 
                 data = results %>% filter( species == 'Y2'), 
                 start = init_pars2, 
                 algorithm = 'port', 
                 lower = lower_HOI_1, 
                 upper = upper_HOI_1)

fit_1_HOI_2 <- nls(form3, 
                  data = results %>% filter( species == 'Y1'), 
                  start = init_pars3, 
                  algorithm = 'port', 
                  lower = lower_HOI_2, 
                  upper = upper_HOI_2)

fit_2_HOI_2 <- nls(form3, 
                  data = results %>% filter( species == 'Y2'), 
                  start = init_pars3, 
                  algorithm = 'port', 
                  lower = lower_HOI_2, 
                  upper = upper_HOI_2)


fit_1_HOI_3 <- nls(form4, 
                   data = results %>% filter( species == 'Y1'), 
                   start = init_pars4, 
                   algorithm = 'port', 
                   lower = lower_HOI_2, 
                   upper = upper_HOI_2)

fit_2_HOI_3 <- nls(form4, 
                   data = results %>% filter( species == 'Y2'), 
                   start = init_pars4, 
                   algorithm = 'port', 
                   lower = lower_HOI_2, 
                   upper = upper_HOI_2)

fit_1
fit_1_HOI_1
fit_1_HOI_2
fit_1_HOI_3

fit_2
fit_2_HOI_1
fit_2_HOI_2
fit_2_HOI_3

results <- 
  results %>% 
  mutate( basic = ifelse( species == 'Y1', predict( fit_1), predict( fit_2))) %>% 
  mutate( HOI_1 = ifelse( species == 'Y1', predict( fit_1_HOI_2), predict( fit_2_HOI_2))) %>% 
  mutate( HOI = ifelse( species == 'Y1', predict( fit_1_HOI_3), predict( fit_2_HOI_3))) 


test <-  
  results %>% 
  gather( fit, y_pred, basic:HOI) 


species_1_fit <- 
  test %>% 
  filter( species == 'Y1') %>%
  filter( fit %in% c('basic', 'HOI')) %>% 
  mutate( lt = as.factor( paste0( B2, fit ))) %>% 
  filter( B2 %in% c(0, 1, 10)) %>% 
  ggplot( aes(x = B1, y = y, color = as.factor(B2))) + 
       geom_point() + 
       geom_line(aes(x = B1, y = y_pred, linetype = fit, group = (lt))) + 
       ylab( 'Fecundity of 1') + 
       xlab( 'Density of 1') + 
  scale_linetype_manual(values = c(2,3)) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  scale_color_discrete('Density of 2')

species_2_fit <- 
  test %>% 
  filter( species == 'Y2') %>%
  filter( fit %in% c('basic', 'HOI')) %>% 
  mutate( lt = as.factor( paste0( B1, fit ))) %>% 
  filter( B1 %in% c(0, 1, 10)) %>% 
  ggplot( aes(x = B2, y = y, color = as.factor(B1))) + 
  geom_point() + 
  geom_line(aes(x = B2, y = y_pred, linetype = fit, group = (lt))) + 
  ylab( 'Fecundity of 2') + 
  xlab( 'Density of 2') + 
  scale_linetype_manual(values = c(2,3)) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  scale_color_discrete('Density of 1')

species_1_fit
species_2_fit

ts_plot

ggsave(ts_plot, file = 'figures/example_timeseries.png', width = 6, height = 5)
ggsave(species_1_fit, file = 'figures/species_1_fit.png', width = 6, height = 5)
ggsave(species_2_fit, file = 'figures/species_2_fit.png', width = 6, height = 5)





