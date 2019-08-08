rm(list = ls())
library(tidyverse)

source('code/plotting_parameters.R')

load('output/model_fits.rda')
source('code/phenomenological_models.R')

pw_comp <- 
  temp_data %>%  
  ungroup () %>% 
  distinct() %>% 
  filter( n_comp < 2 ) %>% 
  select( - starts_with('X'), - time, -h )  %>% 
  gather( comp, density, B1:B3) %>% 
  filter(n_comp == 0 | density > 0) %>% 
  mutate( Competitor = factor(comp, label = species_labs)) %>% 
  mutate( species_lab = factor( species, 
                                labels = c('A) Early', 'B) Mid', 'C) Late')))

pw_grid <- 
  pgrid %>%  
  ungroup () %>% 
  distinct() %>% 
  filter( n_comp < 2 ) %>% 
  gather( comp, density, B1:B3) %>% 
  filter(n_comp == 0 | density > 0) %>% 
  mutate( Competitor = factor(comp, label = species_labs)) %>% 
  mutate( species_lab = factor( species, 
                                labels = c('A) Early', 'B) Mid', 'C) Late')))

pw_grid <- 
  pw_grid %>% 
  gather( model, pred, starts_with('m')) %>% 
  filter( model %in% c('m1', 'm3'))

# plot model predictions and "observed" simulated data ---------------- # 
theme1 <- 
  journal_theme + 
  theme( strip.text = element_text(hjust = 0.1), 
         legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5), 
         legend.key.width = unit(2, 'line'),
         legend.text = element_text(size = 14),
         plot.title = element_text(hjust = 0.5, size = 16))


pw_comp_pred_gg <- 
  pw_comp %>% 
  filter( density < 15) %>% 
  ggplot( aes( x = density, 
               y = y, 
               color = Competitor,
               shape = Competitor)) + 
  geom_point(size = 3) + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3], 
                     name = 'Competitor\nSpecies') +
  scale_shape_manual(values = c(19,17,15), 
                     name = 'Competitor\nSpecies') +  
  scale_x_continuous(breaks = c(0,4,8)) +
  ylab( 'Per Capita Seed Production') + 
  xlab( 'Competitor Density') + 
  ggtitle('Focal Species') +
  guides(color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
         shape = guide_legend(order = 1)) + 
  theme1


pw_comp_pred_gg <- 
  pw_comp_pred_gg +
  geom_line(data = pw_grid, aes( y = pred, linetype = model  ))

pw_comp_pred_gg2 <- 
  pw_comp %>% 
  filter( str_extract( species, pattern = '\\d') != str_extract(comp , pattern = '\\d') ) %>% 
  filter( n_comp == 0 | density > 0) %>% 
  filter( density < 15) %>% 
  ggplot( aes( x = density, y = y, shape = Competitor, color = Competitor)) + 
  geom_point(size = 3) + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3], name = 'Competitor\nSpecies') +
  scale_shape_discrete(name = 'Competitor\nSpecies') + 
  #scale_linetype_discrete(name = 'Best Fit') + 
  scale_x_continuous(breaks = c(0,4,8)) + 
  ylab( 'Per Capita Seed Production') + 
  xlab( 'Competitor Density') + 
  ggtitle('Focal Species')  + 
  guides(shape = guide_legend(order = 1), 
         color = guide_legend(order = 1 )) + #override.aes = list(linetype = 0))) + 
  theme1

temp <- 
  pw_grid %>% 
  filter( model == 'm1') %>%
  filter( str_extract( species, pattern = '\\d') != str_extract(comp , pattern = '\\d')) 


pw_grid %>% 
  filter( comp %in% c('B1', 'B3'), species == 'Y2') %>% 
  filter( model == 'm1')  %>% 
  spread( comp, pred )

pw_comp_pred_gg2 +
  geom_line(data = temp %>% filter( comp == 'B1'), aes( y = pred) ) + 
  geom_line( data = temp %>% filter( comp == 'B3'), aes( y = pred ))
  
pw_grid %>% filter( species == 'Y2', comp == 'B2', model %in% c('m1', 'm3')) %>% spread( model, pred )

ggsave(pw_comp_pred_gg2 +  
         geom_line(data = pw_grid, aes( y = pred, linetype = model  )), 
       filename = 'figures/figure_3.png', 
       width = 7, 
       height = 4)


ggsave(pw_comp_pred_gg, 
       filename = 'figures/figure_S1.png', 
       width = 7, 
       height = 4)

# Compare predicted effects of two species competition to the observed 
# effects of two species competition. Add separate species effects together 
# to predict total effect. Errors are deviations from additivity of 
# competition. 

two_sp_df <- sim_results %>% 
  ungroup() %>%
  filter( n_comp < 3) %>% 
  mutate( species_lab = factor(species, labels = c('A) Early', 'B) Mid', 'C) Late')))

two_sp_df$m2 <- NA
two_sp_df$m2[two_sp_df$species == 'Y1'] <- predict(nls1, newdata = two_sp_df[two_sp_df$species == 'Y1' ,] )
two_sp_df$m2[two_sp_df$species == 'Y2'] <- predict(nls2, newdata = two_sp_df[two_sp_df$species == 'Y2' ,] )
two_sp_df$m2[two_sp_df$species == 'Y3'] <- predict(nls3, newdata = two_sp_df[two_sp_df$species == 'Y3' ,] )

two_sp_df$m1 <- NA
two_sp_df$m1[two_sp_df$species == 'Y1'] <- predict(nls_0_1, newdata = two_sp_df[two_sp_df$species == 'Y1' ,] )
two_sp_df$m1[two_sp_df$species == 'Y2'] <- predict(nls_0_2, newdata = two_sp_df[two_sp_df$species == 'Y2' ,] )
two_sp_df$m1[two_sp_df$species == 'Y3'] <- predict(nls_0_3, newdata = two_sp_df[two_sp_df$species == 'Y3' ,] )

two_sp_df <- 
  two_sp_df %>% 
  gather( mod_type, pred_y, m1:m2 )


theme2 <- 
  journal_theme + 
  theme(strip.text = element_text(hjust = 0.1), 
        legend.title = element_blank(), 
        legend.background = element_blank(), 
        legend.justification = c(1, 1), 
        legend.position = c(1,0.8))


p1 <- 
  two_sp_df %>% 
  filter( mod_type == 'm2') %>% 
  filter( B2 < 15, B3 < 15) %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y1' , B1 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B2, y = y, color = factor(B3), shape = factor(B3) )) + 
  geom_point(size = 3) + 
  geom_line(aes( y = pred_y)) + 
  scale_x_continuous(breaks = c(0,4,8)) + 
  scale_color_manual(values = c('black', 'orange', 'red')) + 
  scale_shape_manual(values = c(1, 0, 2)) + 
  xlab('Mid Competitor\nDensity') + 
  ylab( 'Per Capita Seed Production') + 
  facet_wrap(~ species_lab) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1)) +
  theme2 + 
  annotate(geom = 'text', 
           x = 8, 
           y = 55, 
           label = 'Late Competitor\nDensity', 
           size = 5, 
           hjust = 1)


p2 <- 
  two_sp_df %>% 
  filter( mod_type == 'm2') %>% 
  filter( B1 < 15, B3 < 15) %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y2' , B2 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B3), shape = factor(B3)) ) + 
  geom_point(size = 3) + 
  geom_line(aes( y = pred_y )) + 
  scale_x_continuous(breaks = c(0,4,8)) + 
  scale_color_manual(values = c('black', 'orange', 'red')) + 
  scale_shape_manual(values = c(1, 0, 2)) + 
  xlab('Early Competitor\nDensity') + 
  ylab( 'Per Capita Seed Production') + 
  facet_wrap(~ species_lab) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1)) +
  theme2 + 
  theme(axis.title.y = element_blank()) + 
  annotate(geom = 'text', 
           x = 8, 
           y = 76, 
           label = 'Late Competitor\nDensity', 
           size = 5, 
           hjust = 1)


p3 <- 
  two_sp_df %>% 
  filter( mod_type == 'm2') %>% 
  filter( B1 < 15, B2 < 15) %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y3' , B3 == 0 ) %>% 
  filter( B2 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B2), shape = factor(B2) ) ) + 
  geom_point(size = 3) + 
  geom_line(aes( y = pred_y)) + 
  scale_x_continuous(breaks = c(0,4,8)) + 
  scale_color_manual(values = c('black', 'orange', 'red')) + 
  scale_shape_manual(values = c(1, 0, 2)) + 
  xlab('Early Competitor\nDensity') + 
  ylab( 'Per Capita Seed Production') + 
  facet_wrap(~ species_lab) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1)) +
  theme2 + 
  theme(axis.title.y = element_blank()) + 
  annotate(geom = 'text', 
           x = 8, 
           y = 100, 
           label = 'Mid Competitor\nDensity', 
           size = 5, 
           hjust = 1)


n_comp2 <- grid.arrange(p1, p2, p3, 
                        nrow = 1, 
                        widths = c(0.32, 0.3, 0.3), 
                        top = textGrob('Focal Species', gp = gpar(fontsize = 16)))

ggsave(n_comp2, 
       filename = 'figures/figure_4.png', 
       width = 7.5, 
       height = 4.5)
