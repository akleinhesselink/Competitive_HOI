rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

# graphics themes ------------------------------------------------ # 

journal_theme <- 
  my_theme + 
  theme( axis.title = element_text(size = 12), 
         legend.text = element_text(size = 10), 
         legend.title = element_text(size = 12), 
         strip.text = element_text(size = 12), 
         axis.text = element_text(size = 10))

# Functions ------------------------------------------------------ #

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, r, K, q, m) { B*(q*f(R, r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K) { - sum(B*f(R,r, K)) } # resource (water)

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, r, K, q, m)
    dR <- dRdu(u, B, R, r, K)
    return( list( c(dR, dB))) } )
}

root1 <- function(u, State, parms) with(parms, { State[1] - m[1]*K[1]/(q[1]*r[1]-m[1]) } )
root2 <- function(u, State, parms) with(parms, { State[1] - m[2]*K[2]/(q[2]*r[2]-m[2]) } )
root3 <- function(u, State, parms) with(parms, { State[1] - m[3]*K[3]/(q[3]*r[3]-m[3]) } )

event <- function(u, State, parms) {
  with(parms, {
    State[2:length(State)] <- 0
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

# Run response surface experiments --------------------------- # 

B_init <- expand.grid( 
  B1 = c(0, seq(1, 8, by = 1), 16), 
  B2 = c(0, seq(1, 8, by = 1), 16), 
  B3 = c(0, seq(1, 8, by = 1), 16))

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 )  # filter out three species cases 

out <- list() 
state <- NA
state[1] <- parms$soil_m

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  state[4] <- B_init[i,3]*parms$seedling_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, 200, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root1, 
                  event = list(func = event, root = T), method = 'radau')
  
}

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
results <- data.frame(B_init, results )

pheno <- do.call(rbind, lapply( out, function(x) apply( x, 2, which.max)))
pheno <- data.frame( B_init, pheno )

results <- 
  results %>% 
  mutate( Y1 = X2/(B1*parms$seedling_mass), 
          Y2 = X3/(B2*parms$seedling_mass), 
          Y3 = X4/(B3*parms$seedling_mass)) %>% 
  select( - X1)

results %>% filter( B1 == 1, B2 > 0 , B3 > 0)

results <- 
  results %>% 
  gather( species, y, Y1, Y2, Y3)  %>% 
  mutate( y = y + rnorm(1, 0, 0.001)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, (B1 - 1)*parms$seedling_mass, B1*parms$seedling_mass)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, (B2 - 1)*parms$seedling_mass, B2*parms$seedling_mass)) %>% 
  mutate( B3 = ifelse( str_extract(species, '\\d+') == 3 & B3 > 0, (B3 - 1)*parms$seedling_mass, B3*parms$seedling_mass)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) ) %>% 
  filter( n_comp < 3) %>% 
  filter( !is.na(y))

results %>% filter( B1 == 0, species == 'Y1', B2 > 0 , B3 > 0)
# add HOI column to data ---------------- # 

results <- results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) )

# Filter out NA ------------------------- #  

results <- results %>% 
  filter( !is.na(y))

# Assign lambda -------------------------- # 
results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

# model properties: 

pw_comp_gg <- 
  results %>% 
  ungroup() %>% 
  gather ( comp, density, B1:B3) %>% 
  filter( n_comp < 2 ) %>% 
  filter( density > 0 | y == lambda) %>% 
  ggplot( aes( x = density, y = y, color = comp)) + 
  geom_point() + 
  facet_grid(~species) + 
  scale_color_manual(values = my_colors[1:3]) + 
  ylab( 'Per Capita Growth') + 
  xlab( 'Density') + 
  my_theme

pw_comp_gg + geom_line()

#ggsave(pw_comp_gg + journal_theme, filename = 'ESA_2018/pairwise_comp_no_line.png', width = 10, height = 5.5)

# simple beverton holt with one exponent -------------- # 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3)^tau.'

# species 1 

results %>% head

nls_0_1 <- 
  nls( formula = form_0, 
       data = results %>% 
         filter( species == 'Y1' & n_comp < 2 ), 
       start = list(alpha. = c(2, 2, 2), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 200)

nls_0_1

# species 2 


nls_0_2 <- 
  nls( formula = form_0, 
       data = results  %>% 
         filter( species == 'Y2', n_comp < 2 ) , 
       start = list(alpha. = c(2, 2, 2), tau. = c(2)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 300)

nls_0_2

# species 3 

nls_0_3 <- 
  nls( formula = form_0, 
       data = results %>% 
         filter( species == 'Y3', n_comp < 2 ) , 
       start = list(alpha. = c(2, 2, 2), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 1000)

results$pred_y <- NA
results$pred_y[results$species == 'Y1'] <- predict( nls_0_1, results[ results$species == 'Y1', ])
results$pred_y[results$species == 'Y2'] <- predict( nls_0_2, results[ results$species == 'Y2', ])
results$pred_y[results$species == 'Y3'] <- predict( nls_0_3, results[ results$species == 'Y3', ])


results %>% 
  filter( n_comp < 2 ) %>% 
  gather( comp, density, B1:B3) %>% 
  filter( n_comp == 0 | density > 0) %>% 
  ggplot( aes( x = density, y = y, color = comp)) + 
  geom_point() + 
  geom_line(aes( y = pred_y)) + 
  facet_grid(~species) + 
  scale_color_manual(values = my_colors[1:3]) + 
  ylab( 'Per Capita Growth') + 
  xlab( 'Density') + 
  my_theme

# HOIs in first phase? 
results %>% 
  select ( - starts_with('X')) %>% 
  filter(B1 == 0 , species == 'Y1' ) %>% 
  filter( B3 %in% c(0, 0.005, 0.01, 0.0250 )) %>% 
  ggplot( aes( x = B2, y = y, color = factor(B3) )) + 
  geom_point() + 
  geom_line(aes( y = pred_y) )

results %>% 
  select ( - starts_with('X')) %>% 
  filter(B2 == 0 , species == 'Y2' ) %>% 
  filter( B3 %in% c(0, 0.005, 0.01, 0.0250 )) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B3) )) + 
  geom_point() + 
  geom_line(aes( y = pred_y) )

results %>% 
  select ( - starts_with('X')) %>% 
  filter(B3 == 0 , species == 'Y3' ) %>% 
  filter(B2 %in% c(0, 0.005, 0.01, 0.0250 )) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B2) )) + 
  geom_point() + 
  geom_line(aes( y = pred_y) )


gmean <- function( ... ) prod( ... )^ ( 1/(length( c(...) ) ) )

gmean( 3,2,2,3 ) 
mean(3,2,2,3)                          


B <- rates <- c(1.3,1.2,NA)

B[1] <- 10

for( i in 2:length(rates)) { 
  B[i] <- B[i-1]*(rates[i-1])
}

plot( B )
gmean(rates)
mean(rates)

B[1]*rates[1]*rates[2]
B

