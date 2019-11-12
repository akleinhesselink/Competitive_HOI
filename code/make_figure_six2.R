#
# generate figure 6 in main text -----------------------------------------# 
# 
rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
library(gridExtra)

source('code/plotting_parameters.R')
source('code/simulation_functions.R')

# load parameters --------------------------------- # 
load('output/parms.rda')

# run simulation for mid species effect on late season species---------- # 
R0 <- parms$R0 

seeds_init <- c(0,5,5)
B0 <- seeds_init*parms$seed_mass

state <- c(R0, B0)

out <- ode(state,
           times = seq(0, parms$U, by = 0.001), 
           func = grow, 
           parms = parms, 
           rootfun = root, event = list(func = event, root = T))

tlim <- 85
R <- out[out[,1] < tlim, 2]
g2 <- parms$q*f(R, parms$r[2], parms$K[2]) - parms$m
g3 <- parms$q*f(R, parms$r[3], parms$K[3]) - parms$m
g2[ g2 < 0 ] <- NA
g3[ g3 < 0 ] <- NA

df1 <- data.frame( t = out[out[,1] < tlim, 1], 
                   R = R, 
                   Mid = g2, 
                   Late = g3, 
                   type = 'Early Absent') 

# run simulation with early season species present ---------- # 

seeds_init <- c(5,5,5)
B0 <- seeds_init*parms$seed_mass

state <- c(R0, B0)

out <- ode(state, 
           times = seq(0, 200, by = 0.001), 
           func = grow, 
           parms = parms, 
           rootfun = root, event = list(func = event, root = T))

# organize data for plotting --------------------------------- #

R <- out[out[,1] < tlim, 2]
g2 <- parms$q*f(R, parms$r[2], parms$K[2]) - parms$m
g3 <- parms$q*f(R, parms$r[3], parms$K[3]) - parms$m
g2[ g2 < 0 ] <- NA
g3[ g3 < 0 ] <- NA

df2 <- data.frame( t = out[out[,1] < tlim, 1], 
                   R = R, 
                   Mid = g2, 
                   Late = g3, 
                   type = 'Early Present') 

df <- 
  rbind( df1, df2) %>% 
  group_by( type ) %>% 
  mutate( comp = !is.na(Mid) ) %>% 
  gather( species, rate, c(Mid, Late))  %>% 
  mutate( species = factor(species, levels = c('Mid', 'Late'), ordered = T))

plot_dat1 <- 
  df %>% 
  arrange( species, type, t ) %>% 
  filter( row_number() %% 1000 == 0 | rate < 0.001 ) 

find_spaced <- function(x, n) { 
  
  temp <- 
    data.frame( x = x ) %>%
    mutate( ff = 
              cut( x, breaks = seq( min(x, na.rm = T), max(x, na.rm =T), 
                    length.out = n), 
                    include.lowest = T, 
                    right = T)) %>% 
    group_by( ff ) %>% 
    mutate( keep = x == min(x, na.rm = T)) %>% 
    ungroup() %>% 
    mutate( keep = ifelse( x == max(x, na.rm = T), T, keep))
  
  temp$keep
}


axis_weight <- c(1000, 1)

df_test <- 
  df %>% 
  arrange( species, type, t) %>% 
  group_by( species, type ) %>% 
  mutate( dist = c(0, sqrt( diff(rate*axis_weight[1], 1)^2 + diff(t*axis_weight[2], 1)^2 )), 
          cdist = cumsum(dist)) %>% 
  mutate( plot_points = find_spaced(cdist, n = 30 )) %>% 
  filter( plot_points) %>% 
  mutate(species_comp = paste(species, type, sep = '.')) 

gg1 <- 
  df_test %>%
  ggplot( aes( x = t, y = rate, color = species ) ) + 
  geom_point(aes( shape = species_comp), size = 2) + 
  geom_line(aes( group = species_comp), size = 0.3) + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_shape_manual(values = c(15, 0, 17, 2 )) + 
  xlab( 'Time (d)') + 
  ylab( 'Resource Uptake Rate') + 
  journal_theme + 
  guides( color = F, shape = F) + 
  theme( axis.text.y = element_blank(), 
         legend.position = c(0.25, 0.4), 
         legend.key.width = unit(2, 'line')) + 
  annotate(geom = 'point', 
           c(10, 10), c(0.02, 0.015), size = 3, color = c('red', 'blue'), pch = c(17, 15)) + 
  annotate(geom = 'text', 
           c(20, 20.6), c(0.0205, 0.0151), size = 5, color = c('red', 'blue'), label = c('Mid', 'Late')) + 
  annotate(geom = 'segment', 
           x = c(6, 6), 
           xend = c(14, 14), 
           y = c(0.02, 0.015), 
           yend = c(0.02, 0.015), 
           color = c('red', 'blue')) + 
  annotate(geom = 'text', 
           x = 15, y = 0.028, label = 'Species', size = 6) + 
  ggtitle("A)") + 
  theme(plot.title = element_text(hjust = 0))

activity_bars <- 
  df %>% 
  group_by(species, type ) %>% 
  summarise( xend = max(t[rate > 0], na.rm = T), 
             yend = max(rate, na.rm = T)*0.95) %>%
  mutate( x = xend , y = 0) %>% 
  filter( species == 'Mid')

gg1_notes <- 
  gg1 + 
  annotate( geom = 'segment', 
            activity_bars$x, 
            activity_bars$y, 
            xend = activity_bars$xend, 
            yend = activity_bars$yend, 
            linetype = c(1,2),
            color = 'black', alpha = 0.8) + 
  annotate( geom = 'text', 
            activity_bars$x[1] + 4, 
            activity_bars$yend[1] + 0.003, label = paste('Day',round(activity_bars$x[1])), alpha = 1, size = 4) + 
  annotate( geom = 'text', 
            activity_bars$x[2] - 6, 
            activity_bars$yend[2] + 0.003, label = paste('Day',round(activity_bars$x[2])), alpha = 1, size = 4)



xlims <- ggplot_build(gg1)$layout$panel_scales_x[[1]]$range$range
ylims <- ggplot_build(gg1)$layout$panel_scales_y[[1]]$range$range

gg2 <- 
  df %>% 
  arrange( species, type, t ) %>% 
  filter( comp ) %>% 
  group_by( species, type) %>%
  summarise( avg_rate = mean(rate, na.rm = T)) %>% 
  mutate( type_label = factor( type, label = c('Early\nAbsent','Early\nPresent') )) %>% 
  mutate( species_type = paste( species, type, sep = '.')) %>% 
  ggplot( aes( x = type_label, y = avg_rate, color = species, group = species, shape = species_type)) + 
  geom_point( size = 4) + 
  geom_line() + 
  ylab( 'Avg. Resource Uptake Rate') + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_shape_manual(values = c(15, 0, 17, 2)) + 
  ylim( ylims ) + 
  guides(color = F, shape = F) + 
  journal_theme + 
  theme( axis.text.y = element_blank(), 
         axis.title.x = element_blank(), 
         axis.text.x = element_text(size = 14)) 

gg2 <- 
  gg2 + 
  ggtitle("B)") + 
  theme(plot.title = element_text(hjust = 0))  

gg2 <- 
  gg2 + 
  annotate(geom = 'text', 
           x = c(0.80,0.85), 
           y = c(0.065,0.053), 
           label = c('Mid', 'Late'), 
           color = my_colors[c(2,3)], 
           size = 5) 


gg_both <- 
  grid.arrange(gg1_notes, 
               gg2, 
               nrow = 1, 
               widths = c(0.6, 0.4))

ggsave( 'figures/figure_6.png', 
        gg_both, 
        height = 4, width = 7 )  

