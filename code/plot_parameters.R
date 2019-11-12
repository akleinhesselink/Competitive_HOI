rm(list = ls())

library(tidyverse)
library(gridExtra)
library(grid)

source('code/plotting_parameters.R')
source('code/phenomenological_models.R')
load('output/model_fits.rda')
load('output/processed_results.rda')

theme1 <- 
  journal_theme + 
  theme( legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5), 
         legend.key.width = unit(2, 'line'),
         legend.text = element_text(size = 14),
         plot.title = element_text(hjust = 0.5, size = 20), 
         axis.text.y = element_text(size = 14), 
         axis.text.x = element_blank(), 
         axis.title = element_text(size = 16), 
         plot.margin = unit(c(1,0.5,1,0.5), "lines"), 
         axis.title.x = element_text(margin = margin(2.5,0,0,0,unit = 'lines')), 
         strip.text.y = element_blank(), 
         panel.spacing.y = unit(2, 'lines')) 

model_labs <- c('Hassel', 'Model 2', 'HOI')  


cffs <- lapply( c(fit1, fit2, fit1HOI), coef )

rowids <- expand.grid(Species = species_labs,  
                      Model = model_labs)

cffs_df <- data.frame( rowids, data.frame( do.call( rbind, lapply(cffs, function(x) x[1:5]) ) ))

cffs_df <- 
  cffs_df %>% 
  gather( parameter, value, lambda:tau) %>% 
  bind_rows(
    data.frame( Model = 'HOI', 
            Species = species_labs, 
            do.call( rbind, lapply( cffs, function(x) x[ grep ( 'beta', x = names(x)) ] ) ) ) %>% 
      gather( parameter, value, beta1:beta3))  %>% 
  filter( Model == 'HOI')  %>% 
  filter( !parameter %in% c('lambda', 'tau') ) %>%
  mutate( type = str_extract(parameter, '[a-z]+'), 
          number = str_extract(parameter, '\\d+')) %>%
  mutate( number = factor(number)) 

alpha_labs <- 
  paste0('italic(\u03B1)', 
         paste0( '[', sort( rep( 1:3, 3)), 1:3, ']'))

beta_labs <- 
  paste0( 'italic(\u03B2)', 
          paste0( '[', sort( rep(1:3, 3)), 
                  paste0( '(', c('12', '13', '23'), ')'), ']'))

x_labs <- c(alpha_labs,beta_labs)


label_df <- 
  cffs_df %>% 
  distinct(Species, number, type) %>% 
  arrange(type, Species, number) %>% 
  mutate( label = x_labs) %>% 
  mutate( y_pos = ifelse(type == 'alpha', 
                         -0.1, 
                         -0.1 )) 


letter_df <- 
  label_df %>% 
  distinct(Species, type) %>% 
  mutate( label = paste0 (LETTERS[1:6], ')' ), 
          y_pos = Inf, 
          x_pos = -Inf)


alpha_plot <- 
  cffs_df %>% 
  ggplot( aes( x = number, y = value, color = Species, fill = Species)) + 
  geom_bar(stat = 'identity') + 
  annotate(geom = 'segment', x = 1:3, y = -Inf, xend = 1:3, yend = -0.08 , size = 0.5) + 
  facet_grid(  type ~ Species  ) + 
  scale_fill_manual(name = 'Focal Species', values = my_colors[1:3], guide = F) + 
  scale_color_manual(values = my_colors[1:3], guide = F) +
  ylab( 'Estimate') + 
  xlab( 'Parameter') + 
  ggtitle('Focal Species') +
  theme1 + 
  geom_text(data = label_df, 
            aes( x = number, 
                 y = y_pos, 
                 label = label), 
            parse = T, 
            color = 1, 
            size = 5, 
            vjust = 1) +
  coord_cartesian(ylim = c(0, 1) , clip = 'off') 


alpha_plot <- 
  alpha_plot + 
  geom_text( data = letter_df, 
             aes( x = x_pos, y = y_pos, label = label ), 
             color = 'black', 
             hjust = -0.5, 
             vjust = 1.5 , 
             size = 5) 


alpha_plot %>%   
  ggsave(filename = 'figures/parameter_plot_fig7.png', 
         height = 5, 
         width = 8)

