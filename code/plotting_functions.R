library(tidyverse)
library(deSolve)
library(stringr)
library(gridExtra)
library(scales)

# ggplot themes ------------------------------------------------ # 

my_theme <- 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

journal_theme <- 
  my_theme + 
  theme(  strip.background = element_blank(),
          strip.text = element_text(size = 15),
          axis.text = element_text( size = 12), 
          axis.title = element_text( size = 15), 
          legend.text = element_text( size = 12), 
          legend.title = element_text(size = 15)) 

# colors and labels ----------------------------
my_colors <- c(1,2,4,5)

species_labs <- c('Early', 'Mid', 'Late')

# functions ------------

Rstar <- function(r, K, m, q) { m*K/(q*r-m) }      # resource required for growth to balance loss

plot_timeseries <- function(out, sp_labs = c('Resource', '1','2'), mytheme, fname = 'figures/example_timeseries.png'){ 
  library(gridExtra)
  
  temp <- 
    data.frame(out) %>% 
    gather( var, val, starts_with('X')) %>% 
    mutate( species = var ) %>% 
    mutate( species = factor( species, labels = sp_labs )) %>% 
    filter( time < 150 ) 
    
  resource_plot <- 
    ggplot( temp %>% filter( species == 'Resource'), aes( x = time, y = val )) + 
    geom_line() + 
    ylab( 'Resource') + 
    mytheme + 
    theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.text = element_blank()) 
  
  biomass_plot  <- 
    ggplot( temp %>% filter( species != 'Resource'), aes( x = time, y = val, color = species)) + 
    geom_line() + 
    ylab( 'Biomass') + 
    xlab( 'Day of Year') + 
    scale_color_manual(values = my_colors, guide = F) + 
    mytheme +
    theme(legend.position = c(0.85, 0.5), 
          legend.background = element_rect(colour = 1, size = 0.2), 
          axis.text = element_blank()) 
    

  return(list(resource_plot, biomass_plot))
}

plot_resource_uptake <- function(parms, R = 0:500, spec_labs = c('1','2')){ 
  
  curves <- data.frame(R = R,  mapply(x = as.list(parms$r), y = as.list(parms$K), FUN = function(x, y) { f(R = R, x, y) }) )
  
  curves$X1[curves$R < Rstar(parms$r[1], parms$K[1], parms$m, parms$q) ] <- 0
  curves$X2[curves$R < Rstar(parms$r[2], parms$K[2], parms$m, parms$q) ] <- 0
  curves$X3[curves$R < Rstar(parms$r[3], parms$K[3], parms$m, parms$q) ] <- 0
  
  curves <- 
    curves %>% 
    gather( species, uptake, starts_with('X')) %>% 
    filter( uptake > 0 ) %>% 
    mutate( species = factor(species, labels = spec_labs))
  

  curves %>% 
    ggplot(aes( x = R, y = uptake, color = species )) + 
    geom_line() +
    my_theme +
    ylab( 'Resource uptake rate \n per g per day') + 
    xlab( 'Resource') + 
    scale_color_manual(values = my_colors, 'Species') + 
    theme(legend.position = c(0.75, 0.25), legend.background = element_rect(color = 1, size = 0.25))
  
}

