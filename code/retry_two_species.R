rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, r, K, q, m) { B*(q*f(R, r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K, p, epsilon) { - sum(B*f(R,r, K)) } # resource (water)

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
r <- c(4.2, 2.2) # max uptake rates mm of water per g of plant per day
K <- c(150, 1)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
parms <- list( r = r, K = K, m = m , q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, times = times)

resource_curves <- plot_resource_uptake(parms, spec_labs = c('1', '2'))
resource_curves
ggsave(filename = 'figures/resource_uptake.png', resource_curves, height = 3, width = 4)

R_init <- 200 
seeds_init <- c(1,1)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass)

test <- ode( state, times = 1:200, func = grow, parms = parms)
plot(test)

out <- ode(state, times = 1:200, func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T), method = 'radau')
out

ts_plot <- plot_timeseries(out, sp_labs = c('Resource', '1', '2'))
ts_plot

ggsave( ts_plot, filename = 'figures/example_timeseries.png', height = 4, width = 6)


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
  theme(axis.text = element_blank()) + 
  xlim( 0, 150)

ts_plot

B_init <- expand.grid( B1 = c(0, seq(1, 8, by = 1)), B2 = c(0, seq(1, 8, by = 1)))
B_init <- B_init[-1,]

out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  
  out[[i]] <- ode(state, times = seq(1,200), func = grow, parms = parms, 
                  rootfun = root, event = list(func = event, root = T), method = 'radau')
  
}

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
results
results <- data.frame(B_init, results )

head(results)

results <- 
  results %>% 
  mutate( Y1 = parms$conversion*X2/parms$seedling_mass/B1, Y2 = parms$conversion*X3/parms$seedling_mass/B2) %>% 
  select( - X1)

results <- 
  results %>% 
  gather( species, y, Y1, Y2)  %>% 
  mutate( y = y + rnorm(1, 0, 0.01)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0))

results <- results %>% 
  filter( !is.na(y))

results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

formL <- 'y ~ beta0 + alpha.[1]*B1 + alpha.[2]*B2'
formL_HOI <- 'y ~ beta0 + alpha.[1]*B1 + alpha.[2]*B2 + gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2) + beta.*I(B1*B2)'
formL_HOI_long <- 'y ~ C + alpha.[1]*B1 + alpha.[2]*B2 + gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2) + beta.[1]*I(B1*B2) + beta.[2]*I(B2^3) + beta.[3]*I(B1*B2^2) + beta.[4]*I(B1^2*B2)'


form1 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2)^tau.'
form2 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2))^tau.'
form3 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + beta.*I(B1*B2))^tau.'
form4 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + 
                          gamma.[1]*I(B1^2) + gamma.[2]*I(B2^2) + beta.*I(B1*B2))^tau.'

init_pars1 <- list( alpha. = c(1, 1), tau. = 1)
init_pars2 <- list( alpha. = c(1, 1), gamma. = c(0,0), tau. = 1)
init_pars3 <- list( alpha. = c(1, 1), beta. = 0, tau. = 1)
init_pars4 <- list( alpha. = c(1, 1), gamma. = c(0,0), beta. = 0, tau. = 1)

upper1 <- c(10, 10, 5)
lower1 <- c(0, 0, 0)
upper_HOI_1 <- c(10, 10, 10, 5)
lower_HOI_1 <- c(0, 0, 0, 0, 0)
upper_HOI_2 <- c(100,   100,    100, 100, 100, 5)
lower_HOI_2 <- c(0,       0,      0,   0,   0, 0)

fit_1_L <- nls(formL, 
             data = results %>% filter( species == 'Y1'), 
             start = list(beta0 = 10, alpha. = c(-1,-1)))

fit_2_L <- nls(formL, 
               data = results %>% filter( species == 'Y2'), 
               start = list(beta0 = 10, alpha. = c(-1,-1)))
fit_1_L
fit_2_L

fit_1_L_HOI <- nls(formL_HOI, 
               data = results %>% filter( species == 'Y1'), 
               start = list(beta0 = 10, 
                            alpha. = c(1,1), 
                            gamma. = c(0,0), 
                            beta. = 0))
fit_1_L_HOI

fit_2_L_HOI <- nls(formL_HOI, 
                   data = results %>% filter( species == 'Y2'), 
                   start = list(beta0 = 10, 
                                alpha. = c(1,1), 
                                gamma. = c(0,0), 
                                beta. = 0))
fit_2_L
fit_2_L_HOI

fit_2_L_HOI_long <- nls(formL_HOI_long, 
                   data = results %>% filter( species == 'Y2'), 
                   start = list(C = 10, 
                                alpha. = c(1,1), 
                                gamma. = c(0,0), 
                                beta. = c(0,0,0,0)))

fit_2_L_HOI
fit_2_L_HOI_long




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

fit_2_HOI_3$convInfo

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
  filter( B2 %in% c(0, 1, 8)) %>% 
  ggplot( aes(x = B1, y = y, color = as.factor(B2))) + 
  geom_point() + 
  geom_line(aes(x = B1, y = y_pred, linetype = fit, group = (lt))) + 
  ylab( 'Fecundity of 1') + 
  xlab( 'Density of 1') + 
  scale_linetype_manual(values = c(2,3)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = c(0.8, 0.82)) + 
  scale_color_discrete('Density of 2')


species_2_fit <- 
  test %>% 
  filter( species == 'Y2') %>%
  filter( fit %in% c('basic', 'HOI')) %>% 
  mutate( lt = as.factor( paste0( B1, fit ))) %>% 
  filter( B1 %in% c(0, 1, 8)) %>% 
  ggplot( aes(x = B2, y = y, color = as.factor(B1))) + 
  geom_point() + 
  geom_line(aes(x = B2, y = y_pred, linetype = fit, group = (lt))) + 
  ylab( 'Fecundity of 2') + 
  xlab( 'Density of 2') + 
  scale_linetype_manual(values = c(2,3)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = c(0.8, 0.7)) + 
  guides( color = guide_legend(order = 1), linetype = guide_legend(order = 2)) + 
  scale_color_discrete('Density of 1')

species_1_fit <- 
  species_1_fit + scale_linetype(guide = F)

species_2_fit <- species_2_fit 

both_fits <- grid.arrange(species_1_fit, species_2_fit, ncol = 2)

ggsave(both_fits, file = 'figures/both_fits.png', width = 6, height = 4)

ggsave(species_1_fit, file = 'figures/species_1_fit.png', width = 6, height = 5)
ggsave(species_2_fit, file = 'figures/species_2_fit.png', width = 6, height = 5)


fit_1
fit_1_HOI_1
fit_1_HOI_2
fit_1_HOI_3

fit_2
fit_2_HOI_1
fit_2_HOI_2
fit_2_HOI_3

summary(fit_1)$sigma
test$sigma
summary(fit_2)

fits_df <- data.frame( species = c('1', '2') , 
                       model = c('basic', 'basic'),
            do.call( rbind, lapply( list(fit_1, fit_2), coef) ), 
            residual_error = c( summary(fit_1)$sigma, summary(fit_2)$sigma)) 

fits_df_HOI <- data.frame( species = c('1', '2') , 
                           model = c('HOI', 'HOI'),
                       do.call( rbind, lapply( list(fit_1_HOI_3, fit_2_HOI_3), coef) ), 
                       residual_error = c( summary(fit_1_HOI_3)$sigma, summary(fit_2_HOI_3)$sigma)) 

fit_1
fit_1_HOI_3
fit_2
test <- summary( fit_1_HOI_3)

fits_df <- 
  bind_rows(fits_df, fits_df_HOI ) %>% 
  arrange( species, model ) %>% 
  select( species, model, tau., starts_with('alpha'), starts_with('gamma'), beta., residual_error)


save(fits_df , file = 'data/fits_df.rda')
knitr::kable(fits_df, digits = 3)

fits_df$residual_error[1]/fits_df$residual_error[2]
fits_df$residual_error[3]/fits_df$residual_error[4]

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


# ------------------------- # 

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
which_max <- do.call( rbind, lapply( out, function( x) apply( x, 2, which.max)))

which_max
z <- which_max[, 3]
half_season <- mapply( x = out, y = z, FUN = function(x, y) x[ y, -c(1:2)] )
half_season
half_season <- t( half_season )

out1 <- lapply( out, function(x) x[ 1:70, ] )

results1 <- do.call( rbind, lapply( out1, function(x) apply( x, 2, max )))

results1 <- data.frame(B_init, results1 )
results1

results1 <- 
  results1 %>%
  select( B1, B2, X2:X3) 

results1 <- 
  results1 %>% 
  mutate( Y1 = parms$conversion*X2/parms$seedling_mass/B1, Y2 = parms$conversion*X3/parms$seedling_mass/B2) %>% 
  select( -X2, -X3)

results1 <- 
  results1 %>% 
  gather( species, y, Y1, Y2)  %>% 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0))

results1 <- 
  results1 %>% 
  filter( !is.na(y))

results1 <- 
  results1 %>%
  group_by( species) %>% 
  mutate( lambda = max(y))


fit_1 <- lm(y ~ B1 + B2, 
             data = results1 %>% filter( species == 'Y1'))

summary( fit_1 ) 

fit_2 <- lm(y ~ B1 + B2, 
            data = results1 %>% filter( species == 'Y2'))

summary( fit_2 )



results <- 
  results %>% 
  mutate( pred = ifelse( species == 'Y1', predict(fit_1), NA )) %>% 
  mutate( pred = ifelse( species == 'Y2', predict(fit_2), pred)) 


results %>% head

results %>% 
  filter( species == 'Y1') %>% 
  ggplot( aes( x = B1, y = y, color = as.factor(B2))) + 
  geom_point() + 
  geom_line(aes(x = B1))

results %>% 
  filter( species == 'Y2') %>% 
  ggplot( aes( x = B1, y = y, color = as.factor(B2))) + 
  geom_point() + 
  geom_line(aes(x = B1))


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


summary(fit_2 )
summary(fit_1_HOI_2)
fit_1_HOI_3
fit_2_HOI_3

fit_1
fit_2



