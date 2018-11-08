rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/plotting_functions.R')
source('code/model_functions.R')

# load parameters --------------------------------- # 
load('data/parms.rda')

# run simulation for mid species effect on late season species---------- # 
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

# run simulation with early season species present ---------- # 

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


df %>% mutate( row_number() %% 4)

gg1 <- 
  df %>% 
  arrange( species, type, t ) %>% 
  filter( row_number() %% 10 == 0 ) %>% 
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
  ggtitle("A)") + 
  theme(plot.title = element_text(hjust = 0))

xlims <- ggplot_build(gg1)$layout$panel_scales_x[[1]]$range$range
ylims <- ggplot_build(gg1)$layout$panel_scales_y[[1]]$range$range

gg2 <- 
  df %>% 
  arrange( species, type, t ) %>% 
  filter( row_number() %% 10 == 0 ) %>% 
  filter( comp ) %>% 
  group_by( species, type) %>%
  summarise( avg_rate = mean(rate, na.rm = T)) %>% 
  ggplot( aes( x = type, y = avg_rate, color = species, shape = type)) + 
  geom_point( size = 3) + 
  ylab( 'Avg. Resource Uptake Rate') + 
  xlab( '') + 
  scale_color_manual(values = my_colors[2:3], 'Species') + 
  scale_shape_manual(values = c(19, 1)) + 
  journal_theme + 
  ylim( ylims ) + 
  theme( axis.text.y = element_blank(), axis.title.x = element_text(size = 14)) + 
  guides(color = F, shape = F) 


gg2 <- 
  gg2 + 
  ggtitle("B)") + 
  theme(plot.title = element_text(hjust = 0))
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

gg_both <- grid.arrange(gg1, gg2 + geom_line( aes( group = species)), nrow = 1, widths = c(0.6, 0.4))

ggsave( 'figures/mechanism_of_HOI.png', 
        gg_both, 
        height = 4, width = 7 )  

