rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

# graphics themes ------------------------------------------------ # 

journal_theme <- 
  my_theme + 
  theme( axis.title = element_text(size = 10), 
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
R <- seq(0, 200, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)

R_init <- 200 
seeds_init <- c(0,1,1)
state <- c(R_init, seeds_init)
out <- ode(state, times = seq(0, 200, by = 0.001), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

tlim <- 27
R <- out[out[,1] < tlim, 2]
g2 <- parms$q*f(R, parms$r[2], parms$K[2]) - parms$m
g3 <- parms$q*f(R, parms$r[3], parms$K[3]) - parms$m
g2[ g2 < 0 ] <- NA
g3[ g3 < 0 ] <- NA

df1 <- data.frame( t = out[out[,1] < tlim, 1], R = R, Mid = g2, Late = g3, type = 'Early Absent') 

R_init <- 200 
seeds_init <- c(1,1,1)
state <- c(R_init, seeds_init)
out <- ode(state, times = seq(0, 200, by = 0.001), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

R <- out[out[,1] < tlim, 2]
g2 <- parms$q*f(R, parms$r[2], parms$K[2]) - parms$m
g3 <- parms$q*f(R, parms$r[3], parms$K[3]) - parms$m
g2[ g2 < 0 ] <- NA
g3[ g3 < 0 ] <- NA

df2 <- data.frame( t = out[out[,1] < tlim, 1], R = R, Mid = g2, Late = g3, type = 'Early Present') 

df <- 
  rbind( df1, df2) %>% 
  group_by( type ) %>% 
  mutate( comp = !is.na(Mid) ) %>% 
  gather( species, rate, c(Mid, Late))  %>% 
  mutate( species = factor(species, levels = c('Mid', 'Late'), ordered = T))



gg1 <- 
  df %>% 
  ggplot( aes( x = t, y = rate, color = species, linetype = type)) +
  geom_line() + 
  my_theme + 
  journal_theme + 
  theme( axis.text.y = element_blank(), 
         legend.position = c(0.25, 0.4)) + 
  xlab( 'Time (d)') + 
  ylab( 'Resource Uptake Rate') + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_linetype_manual(values = c(1,2), '')

ylims <- ggplot_build(gg1)$layout$panel_ranges[[1]]$y.range
xlims <- ggplot_build(gg1)$layout$panel_ranges[[1]]$x.range

gg1 <- 
  gg1 + 
  annotate(geom = 'text', label = 'A)', x = xlims[1] + diff(xlims)*0.05 , y = ylims[2] + diff(ylims)*0.05  )

ylims <- ggplot_build(gg1)$layout$panel_ranges[[1]]$y.range
xlims <- ggplot_build(gg1)$layout$panel_ranges[[1]]$x.range

gg2 <- 
  df %>% 
  filter( comp ) %>% 
  group_by( species, type) %>%
  summarise( avg_rate = mean(rate, na.rm = T)) %>% 
  ggplot( aes( x = type, y = avg_rate, color = species, shape = type)) + 
  geom_point( size = 3) + 
  ylab( 'Avg. Resource Uptake Rate') + 
  xlab( '') + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_shape_manual(values = c(19, 1)) + 
  my_theme + 
  journal_theme + 
  ylim( ylims ) + 
  theme( axis.text.y = element_blank(), axis.title.x = element_text(size = 14)) + 
  guides(color = F, shape = F) 


gg2 <- 
  gg2 + 
  annotate(geom = 'text', label = 'B)', x = 0.8, y = ylims[2]) + 
  ylim(ylims)



activity_bars <- 
  df %>% 
  group_by(species, type ) %>% 
  summarise( xend = max(t[rate > 0], na.rm = T), 
             yend = max(rate, na.rm = T)*0.95) %>%
  mutate( x = xend , y = 0) %>% 
  filter( species == 'Mid')


gg1 <- 
  gg1 + 
  annotate( geom = 'segment', 
            activity_bars$x, 
            activity_bars$y, 
            xend = activity_bars$xend, 
            yend = activity_bars$yend, 
            linetype = c(1,2),
            color = 'black', alpha = 0.5) + 
  annotate( geom = 'text', 
            activity_bars$x[1] + 1, 
            activity_bars$yend[1] + 0.003, label = paste('Day',round(activity_bars$x[1])), alpha = 0.5) + 
  annotate( geom = 'text', 
            activity_bars$x[2] - 1, 
            activity_bars$yend[2] + 0.003, label = paste('Day',round(activity_bars$x[2])), alpha = 0.5)


ggsave( 'figures/mechanism_of_HOI.png', 
        grid.arrange(gg1, gg2, nrow = 1, widths = c(0.6, 0.4)), 
        height = 4, width = 6 )  

