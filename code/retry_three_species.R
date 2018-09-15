rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, r, K, q, m) { B*(q*f(R, r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K) { - sum(B*f(R,r, K)) } # resource (water)

dBdudR <- function(R, K, r, q) { K*q*r/(K + R)^2 } # partial derivative of per capita growth with respect to resource concentration
dRdudB <- function(R, K, r){ -r*R/(K + R) } # partial derivative of resource change with respect to biomass

dBdB <- function(R, K1, r1, K2, r2, q) { dBdudR(R, K1, r1, q)*dRdudB(R, K2, r2) } 

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, r, K, q, m)
    dR <- dRdu(u, B, R, r, K)
    return( list( c(dR, dB))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K/(q*r-m) < 0.000001) # logical vector of species to terminate
    State[2:length(State)][ terminate ] <- 0
    return(State)
  })
}

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.6, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(150, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.1         # rate of water evaporation and runoff mm per mm per day
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)

resource_curves <- plot_resource_uptake(parms, spec_labs = c('1', '2', '3'))
resource_curves
ggsave(filename = 'figures/resource_uptake.png', resource_curves, height = 3, width = 4)

R_init <- 200 
seeds_init <- c(1,1, 1)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass, seeds_init[3]*seedling_mass)
test <- ode( state, times = 1:200, func = grow, parms = parms)
plot(test)

out <- ode(state, times = 1:200, func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T), method = 'radau')

ts_plot <- plot_timeseries(out, sp_labs = c('Resource', '1', '2', '3'))
ts_plot
ggsave( ts_plot, filename = 'figures/example_timeseries.png', height = 4, width = 6)

df <- data.frame(out)  
names(df) <- c('time', 'R', '1', '2', '3')

df <- 
  df %>% 
  gather( species, biomass, `1`:`3`) 

ts_plot <- 
  df %>% 
  ggplot( aes( x = time, y = biomass, color = species )) + 
  geom_line() + 
  scale_color_manual(values = my_colors) + 
  my_theme + 
  theme(axis.text = element_blank()) + 
  xlim( 0, 150)

ts_plot

B_init <- expand.grid( 
  B1 = c(0, seq(1, 8, by = 1)), 
  B2 = c(0, seq(1, 8, by = 1)), 
  B3 = c(0, seq(1,8,by=1)))

B_init <- B_init[-1,]

out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  state[4] <- B_init[i,3]*parms$seedling_mass 
  
  out[[i]] <- ode(state, times = seq(1,200), func = grow, parms = parms, 
                  rootfun = root, event = list(func = event, root = T), method = 'radau')
  
}

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
results <- data.frame(B_init, results )
results

pheno <- do.call(rbind, lapply( out, function(x) apply( x, 2, which.max)))
pheno <- data.frame( B_init, pheno )

results <- 
  results %>% 
  mutate( Y1 = parms$conversion*X2/parms$seedling_mass/B1, 
          Y2 = parms$conversion*X3/parms$seedling_mass/B2, 
          Y3 = parms$conversion*X4/parms$seedling_mass/B3) %>% 
  select( - X1)

results <- 
  results %>% 
  gather( species, y, Y1, Y2, Y3)  %>% 
  mutate( y = y + rnorm(1, 0, 0.001)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( B3 = ifelse( str_extract(species, '\\d+') == 3 & B3 > 0, B3 - 1, B3)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) ) %>% 
  filter( n_comp < 3)

# add HOI column to data 

results <- results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) )


results <- results %>% 
  filter( !is.na(y))

results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

formL <- 'y ~ beta0 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3'
formL_HOI <- 'y ~ beta0 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3 + gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2) + gamma.[3]*I(B3^2) + beta.[1]*I(B1*B2) + beta.[2]*I(B1*B3) + beta.[3]*I(B2*B3)'

form1 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3)^tau.'

form2 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3 + 
                      gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2) + gamma.[3]*I(B3^2))^tau.'

form3 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3 + 
                          gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2) + gamma.[3]*I(B3^2) + 
                          beta.[1]*I(B1*B2) + beta.[2]*I(B1*B3) + beta.[3]*I(B2*B3))^tau.'


form_MS_1 <- 'y ~ lambda*exp(-alpha.[1]*B1 + -alpha.[2]*B2 + - alpha.[3]*B3)'

form_1_special <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2] + (alpha.[3]*B3)^tau.[3])'
form_HOI_special <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2] + (alpha.[3]*B3)^tau.[3] + 
                      (beta.[1]*I(B1*B2))^tau.[4] + (beta.[2]*I(B1*B3))^tau.[5] + (beta.[3]*I(B2*B3))^tau.[6])'


init_pars1 <- list( alpha. = c(1, 1, 1), tau. = 1)
init_pars2 <- list( alpha. = c(1, 1, 1), gamma. = c(0,0, 0), tau. = 1)
init_pars3 <- list( alpha. = c(1, 1, 1), gamma. = c(0,0, 0), beta. = c(0,0,0), tau. = 1)

upper1 <- c(10, 10, 10, 5)
lower1 <- c(0, 0, 0, 0)


fit_1_L <- nls(formL, 
               data = results %>% filter( species == 'Y1', n_comp < 3), 
               start = list(beta0 = 10, alpha. = c(-1,-1,-1)))

fit_2_L <- nls(formL, 
               data = results %>% filter( species == 'Y2', n_comp < 3), 
               start = list(beta0 = 10, alpha. = c(-1,-1,-1)))
fit_3_L <- nls(formL, 
               data = results %>% filter( species == 'Y3', n_comp < 3), 
               start = list(beta0 = 10, alpha. = c(-1,-1,-1)))



fit_1_L_HOI <- nls(formL_HOI, 
                   data = results %>% filter( species == 'Y1'), 
                   start = list(beta0 = 10, 
                                alpha. = c(1,1,1), 
                                gamma. = c(0,0,0), 
                                beta. = c(0,0,0)))
fit_2_L_HOI <- nls(formL_HOI, 
                   data = results %>% filter( species == 'Y2'), 
                   start = list(beta0 = 10, 
                                alpha. = c(1,1,1), 
                                gamma. = c(0,0,0), 
                                beta. = c(0,0,0)))
fit_3_L_HOI <- nls(formL_HOI, 
                   data = results %>% filter( species == 'Y3'), 
                   start = list(beta0 = 10, 
                                alpha. = c(1,1,1), 
                                gamma. = c(0,0,0), 
                                beta. = c(0,0,0)))


# fit BH with only one competitor at a time 

small_results <- results %>% filter( (B1 == 0 | B1 > 4),  (B2 == 0 | B2 > 4), (B3 == 0 | B3 > 4) )

fit_3_simple <- nls(y ~ lambda/( 1 + alpha.[1]*B3)^tau., 
                    data = results %>% filter( species == 'Y3', n_comp < 2, B1 == 0 , B2 == 0 ), 
                    start = list(alpha. = 1, tau. = 1))

plot( eval( fit_3_simple$data )$y , predict( fit_3_simple )  )
fit_3_simple

fit_1_simple <- nls(form1, 
                    data = small_results %>% filter( species == 'Y1', n_comp < 2 ), 
                    start = init_pars1, 
                    algorithm = 'port', 
                    lower = lower1, 
                    upper = upper1)

plot( eval( fit_1_simple$data )$y , predict( fit_1_simple )  )

fit_2_simple <- nls(form1, 
                    data = small_results %>% filter( species == 'Y2', n_comp < 2 ), 
                    start = init_pars1, 
                    algorithm = 'port', 
                    lower = lower1, 
                    upper = upper1)

plot( eval( fit_2_simple$data )$y , predict( fit_2_simple )  )


fit_3_simple <- nls(form1, 
                    data = small_results %>% filter( species == 'Y3', n_comp < 2 ), 
                    start = init_pars1, 
                    algorithm = 'port', 
                    lower = lower1, 
                    upper = upper1)

plot( eval( fit_3_simple$data )$y, predict( fit_3_simple )  )

# fit Mayfield Stouffer model with one competitor at a time 

fit_1_MS <- nls(form_MS_1, 
                    data = results %>% filter( species == 'Y1', n_comp < 2 ), 
                    start = list(alpha. = c(0.1, 0.1, 0.1)), 
                    algorithm = 'port', 
                    lower = 0, 
                    upper = 10)
fit_1_MS
plot( eval( fit_1_MS$data )$y, predict( fit_1_MS )  )
abline(0,1)


fit_3_MS <- nls(form_MS_1, 
                data = results %>% filter( species == 'Y3', n_comp < 2 ), 
                start = list(alpha. = c(0.1, 0.1, 0.1)), 
                algorithm = 'port', 
                lower = 0, 
                upper = 10)
fit_3_MS
plot( eval( fit_3_MS$data )$y, predict( fit_3_MS )  )
abline(0,1)



# try Hassel model with varying exponents on each competitive term 
form_1_exp <- 'y ~ lambda/(1 + (alpha.[1]*B2)^tau.[1] + (alpha.[2]*B3)^tau.[2])'
form_1_exp_HOI <- 'y ~ lambda/(1 + (alpha.[1]*B2)^tau.[1] + (alpha.[2]*B3)^tau.[2] + beta.*HOI)'

fit_1_exp <- nls(form_1_exp, 
             data = results %>% filter( species == 'Y1', n_comp < 3, B1 == 0), 
             start = list(alpha. = c(1,1), tau. = c(1,1)), 
             algorithm = 'port', 
             lower = 0, 
             upper = 2)

fit_1_exp_HOI <- nls(form_1_exp_HOI, 
                 data = results %>% filter( species == 'Y1', n_comp < 3, B1 == 0), 
                 start = list(alpha. = c(1,1), tau. = c(1,1), beta. = 1), 
                 algorithm = 'port', 
                 lower = 0, 
                 upper = 2)
fit_1_exp
fit_1_exp_HOI

form_2_exp <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B3)^tau.[2])'
form_2_exp_HOI <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B3)^tau.[2] + beta.*HOI)'

fit_2_exp <- nls(form_2_exp, 
                     data = results %>% filter( species == 'Y2', n_comp < 3, B2 == 0), 
                     start = list(alpha. = c(1,1), tau. = c(1,1)), 
                     algorithm = 'port', 
                     lower = 0, 
                     upper = 2)

fit_2_exp_HOI <- nls(form_2_exp_HOI, 
                 data = results %>% filter( species == 'Y2', n_comp < 3, B2 == 0), 
                 start = list(alpha. = c(1,1), tau. = c(1,1), beta. = 1), 
                 algorithm = 'port', 
                 lower = 0, 
                 upper = 2)

fit_2_exp
fit_2_exp_HOI

form_3_exp <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2])'
form_3_exp_HOI <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2] + beta.*HOI)'

fit_3_exp <- nls(form_3_exp, 
                 data = results %>% filter( species == 'Y3', n_comp < 3, B3 == 0), 
                 start = list(alpha. = c(1,1), tau. = c(1,1)), 
                 algorithm = 'port', 
                 lower = 0, 
                 upper = 2)

fit_3_exp_HOI <- nls(form_3_exp_HOI, 
                 data = results %>% filter( species == 'Y3', n_comp < 3, B3 == 0), 
                 start = list(alpha. = c(1,1), tau. = c(1,1), beta. = -0.5), 
                 algorithm = 'port', 
                 lower = c(0,1,0,0, -2), 
                 upper = c(2,5,1,1,  -1e-6))
fit_3_exp_HOI
fit_3_exp

plot(eval(fit_3_exp_HOI$data)$y, predict(fit_3_exp_HOI))
abline(0,1)

plot( eval(fit_3_exp_HOI$data)$B1*eval(fit_3_exp_HOI$data)$B2 , resid(fit_3_exp_HOI))


form_3_exp_HOI2 <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2] + beta.[1]*HOI + beta.[2]*(I(B1*B2))^0.8 )'

fit_3_exp_HOI2 <- nls(form_3_exp_HOI2, 
                     data = results %>% filter( species == 'Y3', n_comp < 3, B3 == 0), 
                     start = list(alpha. = c(1,1), tau. = c(1,1), beta. = c(-0.5, -0.5)), 
                     algorithm = 'port', 
                     lower = c(0,1,0,0, -40, -2), 
                     upper = c(2,5,1,1, -1e-6, -1e-6))


fit_3_exp
fit_3_exp_HOI
fit_3_exp_HOI2

plot(eval(fit_3_exp_HOI2$data)$B1*eval(fit_3_exp_HOI2$data)$B2,  resid(fit_3_exp_HOI2))



fit_3_df <- data.frame( eval(fit_3_exp$data), basic_mod =  predict( fit_3_exp), HOI_mod = predict(fit_3_exp_HOI))


fit_3_df %>% 
  ggplot( aes( x = B1, y = y, color = factor(B2)) ) + 
  geom_point() + 
  geom_line( aes(y = basic_mod)) + 
  geom_line( aes( y = HOI_mod), linetype = 2)

########## Standard HASSEL model below 

# species 1 as focal
# 2 and 3 as competitors

form1_simple <- 'y ~ lambda/(1 + alpha.[1]*B2 + alpha.[2]*B3)^tau.' 
form1_simple_HOI <- 'y ~ lambda/(1 + alpha.[1]*B2 + alpha.[2]*B3 + beta.[1]*HOI)^tau.'

fit_1 <- nls(form1_simple, 
             data = results %>% filter( species == 'Y1', n_comp < 3, B1 == 0), 
             start = list(alpha. = c(1,1), tau. = 1 ), 
             algorithm = 'port', 
             lower = 0, 
             upper = 2)

fit_1_HOI <- nls(form1_simple_HOI, 
             data = results %>% filter( species == 'Y1', n_comp < 3, B1 == 0), 
             start = list(alpha. = c(1,1), beta. = 1, tau. = 1 ), 
             algorithm = 'port', 
             lower = 0, 
             upper = 2)
fit_1
fit_1_HOI

# species 2 as focal
# 1 and 3 as competitors

form2_simple <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B3)^tau.' 
form2_simple_HOI <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B3 + beta.[1]*HOI)^tau.'

fit_2 <- nls(form2_simple, 
             data = results %>% filter( species == 'Y2', n_comp < 3, B2 == 0), 
             start = list(alpha. = c(1,1), tau. = 1 ), 
             algorithm = 'port', 
             lower = 0, 
             upper = 2)

fit_2_HOI <- nls(form2_simple_HOI, 
                 data = results %>% filter( species == 'Y2', n_comp < 3, B2 == 0), 
                 start = list(alpha. = c(1,1), beta. = 1, tau. = 1 ), 
                 algorithm = 'port', 
                 lower = 0, 
                 upper = 2)
fit_2
fit_2_HOI


# species 3 as focal
# 1 and 2 as competitors


form3_simple <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2)^tau.' 
form3_simple_HOI <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.[1]*HOI)^tau.'

fit_3 <- nls(form3_simple, 
             data = results %>% filter( species == 'Y3', n_comp < 3, B3 == 0), 
             start = list(alpha. = c(1,1), tau. = 1 ), 
             algorithm = 'port', 
             lower = 0, 
             upper = 100)
fit_3

fit_3_HOI <- nls(form3_simple_HOI, 
                 data = results %>% filter( species == 'Y3', n_comp < 3, B3 == 0), 
                 start = list(alpha. = c(1,1), beta. = 1, tau. = 1 ), 
                 algorithm = 'port', 
                 lower = 0, 
                 upper = 100)
fit_3
fit_3_HOI





fit_1 <- nls(form1, 
             data = results %>% filter( species == 'Y1', n_comp < 3), 
             start = init_pars1, 
             algorithm = 'port', 
             lower = lower1, 
             upper = upper1)

fit_2 <- nls(form1, 
             data = results %>% filter( species == 'Y2', n_comp < 3), 
             start = list(alpha. = c(1,1,0), tau. = 1), 
             algorithm = 'port', 
             lower = lower1, 
             upper = upper1)

fit_3 <- nls(form1, 
             data = results %>% filter( species == 'Y3', n_comp < 3), 
             start = list(alpha. = c(0,1,1), tau. = 1), 
             algorithm = 'port', 
             lower = 0, 
             upper = 100)




fit_1_HOI_1 <- nls(form2, 
                   data = results %>% filter( species == 'Y1', n_comp < 3, B1 > 0), 
                   start = init_pars2, 
                   algorithm = 'port', 
                   lower = 0, 
                   upper = 100)

fit_2_HOI_1 <- nls(form2, 
                   data = results %>% filter( species == 'Y2', n_comp < 3), 
                   start = init_pars2, 
                   algorithm = 'port', 
                   lower = 0, 
                   upper = 100)

fit_3_HOI_1 <- nls(form2, 
                   data = results %>% filter( species == 'Y3', n_comp < 3), 
                   start = init_pars2, 
                   algorithm = 'port', 
                   lower = 0, 
                   upper = 100)


fit_1_HOI_2 <- nls(form3, 
                   data = results %>% filter( species == 'Y1', n_comp < 3), 
                   start = init_pars3, 
                   algorithm = 'port', 
                   lower = 0, 
                   upper = 100)

coefs <- coef(fit_1_HOI_2)
alphas <- grep ('alpha.',  names( coefs))
gammas <- grep ('gamma.',  names( coefs))
betas  <- grep ('beta.', names( coefs))
tau    <- grep ('tau.', names(coefs))

init_pars <- list( alpha. = coefs[alphas], 
                   gamma. = coefs[gammas], 
                   beta. = coefs[betas], 
                   tau. = coefs[tau] )

fit_1_HOI_2 <- nls(form3, 
                   data = results %>% filter( species == 'Y1', n_comp < 3), 
                   start = init_pars, 
                   algorithm = 'port', 
                   lower = c(rep(0,3), c(-1,-1,0), c(-1,0,0), 0), 
                   upper = 100)


fit_2_HOI_2 <- nls(form3, 
                   data = results %>% filter( species == 'Y2', n_comp < 3), 
                   start = init_pars3, 
                   algorithm = 'port', 
                   lower = c(0, 0), 
                   upper = 100)

coefs <- coef(fit_2_HOI_2)
alphas <- grep ('alpha.',  names( coefs))
gammas <- grep ('gamma.',  names( coefs))
betas  <- grep ('beta.', names( coefs))
tau    <- grep ('tau.', names(coefs))

init_pars <- list( alpha. = coefs[alphas], 
                   gamma. = coefs[gammas], 
                   beta. = coefs[betas], 
                   tau. = coefs[tau] )

fit_2_HOI_2 <- nls(form3, 
                   data = results %>% filter( species == 'Y2', n_comp < 3), 
                   start = init_pars3, 
                   algorithm = 'port', 
                   lower = c(rep(0,3), c(-1,0,0), c(0,0,0), 0), 
                   upper = 100)


fit_3_HOI_2 <- nls(form3, 
                   data = results %>% filter( species == 'Y3', n_comp < 3 ), 
                   start = init_pars3, 
                   algorithm = 'port', 
                   lower = c(rep(0,2), -10,-1,-1,rep(0,4), 0), 
                   upper = 1000)

fit_3_HOI_2

coefs <- coef(fit_3_HOI_2)
alphas <- grep ('alpha.',  names( coefs))
gammas <- grep ('gamma.',  names( coefs))
betas  <- grep ('beta.', names( coefs))
tau    <- grep ('tau.', names(coefs))

init_pars <- list( alpha. = coefs[alphas],
                   gamma. = coefs[gammas],
                   beta. = coefs[betas], 
                   tau. = coefs[tau] )

fit_3_HOI_2 <- nls(form3, 
                   data = results %>% filter( species == 'Y3', n_comp < 3), 
                   start = init_pars, 
                   algorithm = 'port', 
                   lower = c(rep(0,2), -10,-1,-1,rep(0,4), 0), 
                   upper = 1000)


fit_1
fit_1_HOI_1
fit_1_HOI_2

fit_2
fit_2_HOI_1
fit_2_HOI_2

fit_3
fit_3_HOI_1
fit_3_HOI_2

results_list <- split( results, results$species)
predict_fit <- mapply(list(fit_1, fit_2, fit_3), results_list, FUN = predict)
results_list <- mapply( results_list, as.list(data.frame(predict_fit)), FUN = function(x, y) {x$basic <- y; return(x) } , SIMPLIFY = F)
plot( results_list$Y1$y, results_list$Y1$basic)
plot( results_list$Y2$y, results_list$Y2$basic)
plot( results_list$Y3$y, results_list$Y3$basic)

results <- do.call(rbind, results_list)

results_list <- split( results, results$species)
predict_fit <- mapply(list(fit_1_HOI_2, fit_2_HOI_2, fit_3_HOI_2), results_list, FUN = predict)
results_list <- mapply( results_list, as.list(data.frame(predict_fit)), FUN = function(x, y) {x$HOI <- y; return(x) } , SIMPLIFY = F)
plot( results_list$Y1$y, results_list$Y1$HOI)
plot( results_list$Y2$y, results_list$Y2$HOI)
plot( results_list$Y3$y, results_list$Y3$HOI)

results <- do.call(rbind, results_list)

test <- 
  results %>% 
  ungroup() %>% 
  gather( model, y_pred, basic, HOI) %>% 
  select( species, B1:B3, y, model, y_pred)

my_theme1 <- my_theme + theme(legend.position = c(0.5, 0.9), legend.justification = c(0,1)) 

gg_112 <- 
  test %>% 
  filter( species == 'Y1', B2 %in% c(0, 1, 8), B3 == 0) %>% 
  ggplot( aes( x = B1, y = y, color = as.factor( B2) )) + 
  geom_point() + 
  geom_line( aes( x = B1, y = y_pred, linetype = model)) + 
  ylab( 'Fecundity of 1') + 
  xlab( 'Density of 1') + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  scale_color_discrete('Density of 2') + 
  my_theme1

gg_113 <- 
  test %>% 
  filter( species == 'Y1', B3 %in% c(0, 1, 8), B2 == 0) %>% 
  ggplot( aes( x = B1, y = y, color = as.factor( B3) )) + 
  geom_point() + 
  geom_line( aes( x = B1, y = y_pred, linetype = model)) + 
  ylab( '') + 
  xlab( 'Density of 1') + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  scale_color_discrete('Density of 3') + 
  my_theme1

gg_123 <- 
  test %>% 
  filter( species == 'Y1', B3 %in% c(0, 1, 8), B1 == 0) %>% 
  ggplot( aes( x = B2, y = y, color = as.factor( B3) )) + 
  geom_point() + 
  geom_line( aes( x = B2, y = y_pred, linetype = model)) + 
  ylab( '') + 
  xlab( 'Density of 2') + 
  scale_linetype_manual(values = c(2,3)) + 
  guides( color = guide_legend(order = 1), linetype = guide_legend(order = 2)) + 
  scale_color_discrete('Density of 3') + 
  my_theme1

species_1_fits <- grid.arrange(gg_112, gg_113, gg_123, ncol = 3)

gg_221 <- 
  test %>% 
  filter( species == 'Y2', B1 %in% c(0, 1, 8), B3 == 0) %>% 
  ggplot( aes( x = B2, y = y, color = as.factor( B1) )) + 
  geom_point() + 
  geom_line( aes( x = B2, y = y_pred, linetype = model)) + 
  ylab( 'Fecundity of 2') + 
  xlab( 'Density of 2') + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  scale_color_discrete('Density of 1') + 
  my_theme1


gg_223 <- 
  test %>% 
  filter( species == 'Y2', B3 %in% c(0, 1, 8), B1 == 0) %>% 
  ggplot( aes( x = B2, y = y, color = as.factor( B3) )) + 
  geom_point() + 
  geom_line( aes( x = B2, y = y_pred, linetype = model)) + 
  ylab( '') + 
  xlab( 'Density of 2') + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  scale_color_discrete('Density of 3') + 
  my_theme1

gg_213 <- 
  test %>% 
  filter( species == 'Y2', B3 %in% c(0, 1, 8), B2 == 0) %>% 
  ggplot( aes( x = B1, y = y, color = as.factor( B3) )) + 
  geom_point() + 
  geom_line( aes( x = B1, y = y_pred, linetype = model)) + 
  ylab( '') + 
  xlab( 'Density of 1') + 
  scale_linetype_manual(values = c(2,3)) + 
  guides( color = guide_legend(order = 1), linetype = guide_legend(order = 2)) + 
  scale_color_discrete('Density of 3') + 
  my_theme1

species_2_fits <- grid.arrange(gg_221, gg_223, gg_213, ncol = 3)

gg_331 <- 
  test %>% 
  filter( species == 'Y3', B1 %in% c(0, 1, 8), B2 == 0) %>% 
  ggplot( aes( x = B3, y = y, color = as.factor( B1) )) + 
  geom_point() + 
  geom_line( aes( x = B3, y = y_pred, linetype = model)) + 
  ylab( 'Fecundity of 3') + 
  xlab( 'Density of 3') + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  scale_color_discrete('Density of 1') + 
  my_theme1

gg_332 <- 
  test %>% 
  filter( species == 'Y3', B2 %in% c(0, 1, 8), B1 == 0) %>% 
  ggplot( aes( x = B3, y = y, color = as.factor(B2) )) + 
  geom_point() + 
  geom_line( aes( x = B3, y = y_pred, linetype = model)) + 
  ylab( '') + 
  xlab( 'Density of 3') + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  scale_color_discrete('Density of 2') + 
  my_theme1

gg_312 <- 
  test %>% 
  filter( species == 'Y3', B2 %in% c(0, 1, 8), B3 == 0) %>% 
  ggplot( aes( x = B1, y = y, color = as.factor( B2) )) + 
  geom_point() + 
  geom_line( aes( x = B1, y = y_pred, linetype = model)) + 
  ylab( '') + 
  xlab( 'Density of 1') + 
  scale_linetype_manual(values = c(2,3)) + 
  guides( color = guide_legend(order = 1), linetype = guide_legend(order = 2)) + 
  scale_color_discrete('Density of 2') + 
  my_theme1


species_3_fits <- grid.arrange(gg_331, gg_332, gg_312, ncol = 3)

gg_123 + ylab('Fecundity of 1') + guides( linetype = F) + annotate(geom = 'text', label = 'a', x = 0, y = 70)

ggsave( grid.arrange( gg_123 + ylab('Fecundity of 1') + 
                        guides( linetype = F) + 
                        annotate(geom = 'text', label = 'a', x = 0, y = 70), 
              gg_213 + 
                ylab('Fecundity of 2') + 
                guides( linetype = F) + 
                annotate(geom = 'text', label = 'b', x = 0, y = 90), 
              gg_312 + ylab('Fecundity of 3') + 
                annotate(geom = 'text', label = 'c', x = 0, y = 110), 
              ncol = 3, nrow = 1), 
        file = 'figures/threeway_fits.png', width = 7.5, height = 4)

ggsave(species_1_fits, file = 'figures/species_1_fits_three_species.png', width = 7.5, height = 4)
ggsave(species_2_fits, file = 'figures/species_2_fits_three_species.png', width = 7.5, height = 4)
ggsave(species_3_fits, file = 'figures/species_3_fits_three_species.png', width = 7.5, height = 4)

saveRDS(test, 'data/simulation_data.rds')




fits_df <- data.frame( species = c('1', '2', '3') , 
                       model = 'basic',
                       do.call( rbind, lapply( list(fit_1, fit_2, fit_3), coef) ), 
                       residual_error = c( summary(fit_1)$sigma, summary(fit_2)$sigma, summary(fit_3)$sigma) ) 

fits_df_HOI <- data.frame( species = c('1', '2', '3') , 
                           model = 'HOI',
                           do.call( rbind, lapply( list(fit_1_HOI_2, fit_2_HOI_2, fit_3_HOI_2), coef) ), 
                           residual_error = c( summary(fit_1_HOI_2)$sigma, summary(fit_2_HOI_2)$sigma, summary(fit_3_HOI_2)$sigma)) 


fits_df <- 
  bind_rows(fits_df, fits_df_HOI ) %>% 
  arrange( species, model ) %>% 
  select( species, model, tau., starts_with('alpha'), starts_with('gamma'), starts_with('beta.'), residual_error)

fits_df
save(fits_df , file = 'data/fits_df.rda')
knitr::kable(fits_df, digits = 2)
