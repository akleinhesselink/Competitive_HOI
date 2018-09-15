rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

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
results <- results %>% mutate( exp_id = row_number())

pheno <- do.call(rbind, lapply( out, function(x) apply( x, 2, which.max)))
pheno <- data.frame( B_init, pheno )
nrow(results)
nrow(pheno)
length(out)

dRdu_vec <- Vectorize(FUN = dRdu)

sp2_R <- sp3_R <- sp1_R <- list()

for( i in 1:length(out)){ 
  sp1_R[[i]] <- dRdu_vec( u = out[[i]][,1], B = out[[i]][, 3], R = out[[i]][,2], r = parms$r[1], K = parms$K[1])
  sp2_R[[i]] <- dRdu_vec( u = out[[i]][,1], B = out[[i]][, 4], R = out[[i]][,2], r = parms$r[2], K = parms$K[2])
  sp3_R[[i]] <- dRdu_vec( u = out[[i]][,1], B = out[[i]][, 5], R = out[[i]][,2], r = parms$r[3], K = parms$K[3])
}

sp_R_U <- mapply( sp1_R, sp2_R, sp3_R, FUN = function(Y1,Y2,Y3) { data.frame( Y1, Y2, Y3)}, SIMPLIFY = F)

phase_df <- pheno[ , c('X2', 'X3', 'X4')] 

R_stars <- Rstar(r = parms$r, K = parms$K, m = parms$m, q = parms$q)

sp1_R[[1]]
out[[1]][, 2] > R_stars[1]

for( i in 1:nrow(phase_df)){ 
  sp_R_U[[i]]$phase1 <- 1:nrow( sp_R_U[[i]]) > 1 & out[[i]][,2] > R_stars[1] 
  sp_R_U[[i]]$phase2 <- 1:nrow( sp_R_U[[i]]) > 1 & !sp_R_U[[i]]$phase1 & out[[i]][,2] > R_stars[2]
  sp_R_U[[i]]$phase3 <- 1:nrow( sp_R_U[[i]]) > 1 & !sp_R_U[[i]]$phase1 & !sp_R_U[[i]]$phase2 & out[[i]][,2] > R_stars[3]
}

R_use_summary <- list()

for( i in 1:length(sp_R_U)){ 
  R_use_summary[[i]] <- 
    sp_R_U[[i]] %>% 
    gather( phase, val, phase1:phase3) %>% 
    filter( val ) %>% 
    select( -val ) %>% 
    gather( species, dRdu, Y1:Y3) %>% 
    group_by( phase, species ) %>%
    summarise( R_phase = sum( dRdu) ) %>% 
    arrange( species, phase) %>% 
    group_by( species) %>% 
    mutate( R_total = sum(R_phase), prop_phase = R_phase/R_total) %>% 
    select(-R_phase, -R_total) %>% 
    spread(phase, prop_phase, fill = 0) %>% 
    mutate( exp_id = i)
}

R_use_summary <- do.call(bind_rows, R_use_summary)
R_use_summary <- R_use_summary %>% select( exp_id, species, phase1, phase2, phase3)


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
  filter( n_comp < 3 )

# add HOI column to data 
results <- results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) )

results <- results %>% 
  filter( !is.na(y))

results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

results %>% head


test <- left_join(results, R_use_summary, by = c('exp_id', 'species'))
test %>% head

test$species

test %>% 
  select( -starts_with('X'), -time) %>% 
  filter( species == 'Y3', B3 == 0, B1 %in% c(0,3), B2 %in% c(0,3)) %>% 
  gather( phase, R, phase1:phase3) %>% 
  mutate( B2_pres = B2 > 0 , B1_pres = B1 > 0) %>% 
  ggplot( aes( x = phase, y = R)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  facet_grid(B2_pres ~ B1_pres  )


test %>% 
  select( -starts_with('X'), -time) %>% 
  filter( species == 'Y2', B1 %in% c(0,3), B2 == 0, B3 == 0) %>% 
  gather( phase, R, phase1:phase3) %>% 
  group_by( B1 == 3, phase) %>% 
  summarise( mean(R)) %>% 
  ggplot( aes( x = phase, y = `mean(R)`, fill = `B1 == 3`)) + geom_bar(stat = 'identity', position = 'dodge')

