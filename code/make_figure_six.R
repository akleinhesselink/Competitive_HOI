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

gg1 <- 
  df %>% 
  arrange( species, type, t ) %>% 
  filter( row_number() %% 100 == 0 | rate < 0.0001 ) %>% 
  ggplot( aes( x = t, y = rate, color = species, linetype = type)) +
  geom_line() + 
  xlab( 'Time (d)') + 
  ylab( 'Resource Uptake Rate') + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_linetype_manual(values = c(1,2), '') + 
  journal_theme + 
  theme( axis.text.y = element_blank(), 
         legend.position = c(0.25, 0.4), 
         legend.key.width = unit(2, 'line'))

gg1 <- 
  gg1 + 
  ggtitle("A)") + 
  theme(plot.title = element_text(hjust = 0))

xlims <- ggplot_build(gg1)$layout$panel_scales_x[[1]]$range$range
ylims <- ggplot_build(gg1)$layout$panel_scales_y[[1]]$range$range

gg2 <- 
  df %>% 
  arrange( species, type, t ) %>% 
  filter( row_number() %% 100 == 0 ) %>% 
  filter( comp ) %>% 
  group_by( species, type) %>%
  summarise( avg_rate = mean(rate, na.rm = T)) %>% 
  mutate( type_label = factor( type, label = c('Early\nAbsent','Early\nPresent') )) %>% 
  ggplot( aes( x = type_label, y = avg_rate, color = species, shape = type_label)) + 
  geom_point( size = 3) + 
  ylab( 'Avg. Resource Uptake Rate') + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_shape_manual(values = c(19, 1)) + 
  ylim( ylims ) + 
  guides(color = F, shape = F) + 
  journal_theme + 
  theme( axis.text.y = element_blank(), 
         axis.title.x = element_blank(), 
         axis.text.x = element_text(size = 14)) 

gg2 <- 
  gg2 + 
  geom_line( aes( group = species)) + 
  ggtitle("B)") + 
  theme(plot.title = element_text(hjust = 0))  

gg2 <- 
  gg2 + 
  annotate(geom = 'text', 
           x = c(0.9,0.9), 
           y = c(0.065,0.052), 
           label = c('Mid', 'Late'), 
           color = my_colors[c(2,3)], 
           size = 5) 

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
            activity_bars$x[1] + 5, 
            activity_bars$yend[1] + 0.003, label = paste('Day',round(activity_bars$x[1])), alpha = 1, size = 4) + 
  annotate( geom = 'text', 
            activity_bars$x[2] - 5, 
            activity_bars$yend[2] + 0.003, label = paste('Day',round(activity_bars$x[2])), alpha = 1, size = 4)

gg_both <- 
  grid.arrange(gg1_notes, 
               gg2, 
               nrow = 1, 
               widths = c(0.6, 0.4))

ggsave( 'figures/figure_6.png', 
        gg_both, 
        height = 4, width = 7 )  

