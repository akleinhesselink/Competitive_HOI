rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/plotting_functions.R')

load( 'data/sim_results.rda')
load( 'data/parms.rda')

sim_results <- 
  sim_results %>% 
  mutate( Y1 = parms$conversion*X2/parms$seedling_mass/B1, 
          Y2 = parms$conversion*X3/parms$seedling_mass/B2, 
          Y3 = parms$conversion*X4/parms$seedling_mass/B3) %>% 
  select( - X1)

sim_results <- 
  sim_results %>% 
  gather( species, y, Y1, Y2, Y3)  %>% 
  mutate( y = y + rnorm(1, 0, 0.001)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( B3 = ifelse( str_extract(species, '\\d+') == 3 & B3 > 0, B3 - 1, B3)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) ) %>% 
  filter( n_comp < 3)

# add HOI column to data ---------------- # 

sim_results <- sim_results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) )


# Filter out NA ------------------------- #  

sim_results <- sim_results %>% 
  filter( !is.na(y))

# Assign lambda -------------------------- # 
sim_results <- sim_results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

# Organize simulation results------------- #  

pw_comp_df <- 
  sim_results %>% 
  ungroup() %>%
  filter( n_comp < 2) %>% 
  gather( comp, density, B1:B3) %>% 
  filter( n_comp == 0 | density > 0) %>% 
  mutate( Competitor = factor(comp, label = species_labs)) %>% 
  mutate( species_lab = factor(species, labels = species_labs ) ) 

pw_comp_gg <- 
  pw_comp_df %>% 
  filter( species_lab != Competitor) %>%
  filter( density < 15) %>% 
  ggplot( aes( x = density, y = y, color = Competitor)) + 
  geom_point() + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3]) + 
  ylab( 'Per Capita Fecundity') + 
  xlab( 'Density') + 
  my_theme

# Fit Simple Beverton-Holt/Hassel Model -------------- # 

# fit species 1  ------------------ # 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3)^tau.'

pw_comp_df1 <- pw_comp_df %>% 
  filter( species_lab != Competitor) %>% 
  filter( species == 'Y1') %>% 
  spread(comp, density, fill = 0)

form_0 <- 'y ~ lambda/(1 + alpha.[1]*B2 + alpha.[2]*B3)^tau.'

nls_0_1 <- 
  nls( formula = form_0, 
       data = pw_comp_df1 , 
       start = list(alpha. = c(0.7,0.3), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

pw_comp_df1$pred_0 <- predict( nls_0_1 )

# fit species 2  ------------------------- # 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B3)^tau.'

pw_comp_df2 <- pw_comp_df %>% 
  filter( species_lab != Competitor) %>% 
  filter( species == 'Y2') %>% 
  spread(comp, density, fill = 0)

nls_0_2 <- 
  nls( formula = form_0, 
       data = pw_comp_df2 , 
       start = list(alpha. = c(1,0.3), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

pw_comp_df2$pred_0 <- predict( nls_0_2 ) 

# Fit species 3 ---------------------------- # 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2)^tau.'

pw_comp_df3 <- pw_comp_df %>% 
  filter( species_lab != Competitor) %>% 
  filter( species == 'Y3') %>% 
  spread(comp, density, fill = 0)

nls_0_3 <- 
  nls( formula = form_0, 
       data = pw_comp_df3 , 
       start = list(alpha. = c(1, 0.7), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 50)

pw_comp_df3$pred_0 <- predict( nls_0_3)

# Fit second model, modified Beverton-Holt/Hassel. ----------- # 
# Competitor densities contribute non-linearly to competition. 

# fit species 1 --------------------- # 

form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B2)^tau.[1] + (alpha.[2]*B3)^tau.[2] )'

nls1 <- 
  nls( formula = form_1, 
       data = pw_comp_df1 , 
       start = list(alpha. = c(0.7, 0.3), tau. = c( 1,1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

pw_comp_df1$pred_1 <- predict( nls1 )

# species 2 ------------------------- # 

form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B3)^tau.[2] )'

nls2 <- 
  nls( formula = form_1, 
       data = pw_comp_df2 , 
       start = list(alpha. = c(1,0.3), tau. = c(1,1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

pw_comp_df2$pred_1 <- predict( nls2 )

# species 3 ------------------------- # 
form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2] )'

nls3 <- 
  nls( formula = form_1, 
       data = pw_comp_df3 , 
       start = list(alpha. = c(1, 0.7), tau. = c(1, 1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 5)

pw_comp_df3$pred_1 <- predict( nls3 )

# Plot model fits for each species   --------------------------- # 

pw_comp_df <- do.call(rbind, list( pw_comp_df1 %>% gather( comp, density, starts_with('B')), 
                                   pw_comp_df2 %>% gather( comp, density, starts_with('B')), 
                                   pw_comp_df3 %>% gather( comp, density, starts_with('B'))))

pw_comp_df <- 
  pw_comp_df %>% 
  gather( pred_type, pred, pred_0:pred_1)

pw_comp_df <- 
  pw_comp_df %>% 
  mutate( Model = factor( pred_type, labels = c('1', '2') )) %>%
  mutate( species_lab = factor( species_lab, labels = c('A) Early', 'B) Mid', 'C) Late')))

pw_comp_pred_gg <- 
  pw_comp_df %>% 
  mutate( `Best Fit` = factor(Model, labels = c('model 1', 'model 2'))) %>% 
  filter( n_comp == 0 | density > 0) %>% 
  filter( density < 15) %>% 
  ggplot( aes( x = density, y = y, 
               color = Competitor,
               shape = Competitor)) + 
  geom_point(size = 3) + 
  geom_line(aes( y = pred, linetype = `Best Fit` )) + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3], 
                     name = 'Competitor') +
  scale_shape_manual(values = c(19,17,15), 
                     name = 'Competitor') +  
  scale_linetype_discrete(name = 'Best Fit') + 
  ylab( 'Per Capita Fecundity') + 
  xlab( 'Competitor Density') + 
  my_theme + 
  journal_theme + 
  guides(color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
         shape = guide_legend(order = 1, override.aes = list(linetype = 0)), 
         linetype = guide_legend(order = 2)) + 
  theme( strip.text = element_text(hjust = 0.1), 
         legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5), 
         legend.key.width = unit(2, 'line'))

pw_comp_pred_gg2 <- 
  pw_comp_df %>% 
  mutate( `Best Fit` = factor(Model, labels = c('Model 1', 'Model 2'))) %>% 
  filter( n_comp == 0 | density > 0) %>% 
  filter( density < 15, Model == 2) %>% 
  ggplot( aes( x = density, y = y, shape = Competitor, color = Competitor)) + 
  geom_point(size = 3) + 
  geom_line(aes( y = pred, linetype = `Best Fit` )) + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3], name = 'Competitor') +
  scale_shape_discrete(name = 'Competitor') + 
  scale_linetype_discrete(name = 'Best Fit') + 
  ylab( 'Per Capita Fecundity') + 
  xlab( 'Competitor Density') + 
  my_theme + 
  journal_theme + 
  guides(shape = guide_legend(order = 1), 
         color = guide_legend(order = 1, 
                              override.aes = list(linetype = 0))) + 
  theme( strip.text = element_text(hjust = 0.1), 
         legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5), 
         legend.key.width = unit(2, 'line'))

pw_comp_pred_gg <- 
  pw_comp_pred_gg + 
  scale_x_continuous(breaks = c(0,4,8))

pw_comp_pred_gg2 <- 
  pw_comp_pred_gg2 + 
  scale_x_continuous(breaks = c(0,4,8)) 

ggsave(pw_comp_pred_gg, filename = 'figures/appendix_pairwise_comp_with_line.png', width = 7, height = 4)

ggsave(pw_comp_pred_gg2, filename = 'figures/pairwise_comp_with_line.png', width = 7, height = 4)

# Compare predicted effects of two species competition to the observed effects of two species competition. 
# Use the model fits above and add species effects together--i.e. deviation from additivity. 

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
  xlab('Density of\nMid Species') + 
  ylab( 'Per Capita Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme( legend.title = element_blank(), 
         legend.position = c(0.7, 0.8), 
        strip.text = element_text(hjust = 0.1)) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1, override.aes = list(linetype = 0))) +
  annotate(geom = 'text', x = 8, y = 55, label = 'Density of\nLate Species', size = 5, hjust = 1)


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
  xlab('Density of\nEarly Species') + 
  ylab( 'Fecundity') + 
  facet_wrap(~ species_lab) + 
  journal_theme + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.7, 0.8), 
        strip.text = element_text(hjust = 0.1)) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1, override.aes = list(linetype = 0))) +
  annotate(geom = 'text', x = 8, y = 76, label = 'Density of\nLate Species', size = 5, hjust = 1)


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
  xlab('Density of\nEarly Species') + 
  ylab( 'Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme( legend.title = element_blank()) + 
  theme(legend.position = c(0.7, 0.8), 
        strip.text = element_text(hjust = 0.1)) + 
  guides( 
    color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
    shape = guide_legend(order = 1, override.aes = list(linetype = 0))) +
  annotate(geom = 'text', x = 8, y = 100, label = 'Density of\nMid Species', size = 5, hjust = 1)


temp_legend_theme <- 
  theme(legend.background = element_blank(), 
      legend.title = element_blank(), 
      legend.justification = c(1, 1), 
      legend.position = c(1,0.8))

p1 <- 
  p1 + 
  temp_legend_theme

p2 <- 
  p2 + 
  temp_legend_theme + 
  theme(axis.title.y = element_blank())

p3 <- 
  p3 + 
  temp_legend_theme + 
  theme(axis.title.y = element_blank())


n_comp2 <- grid.arrange(p1, p2, p3, 
                        nrow = 1, 
                        widths = c(0.32, 0.3, 0.3))

ggsave(n_comp2, 
       filename = 'figures/two_sp_comp_pw_line.png', 
       width = 7, 
       height = 4)




# Plot average error in two species communities ---------------------------------------# 
RMSE_plot_both_models <- 
  two_sp_df %>% 
  select( -time , -c(X2:X4), -lambda) %>%
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>% 
  filter( n_comp > 0 ) %>% 
  group_by(species, mod_type, HOI) %>%
  summarise( RMSE = sqrt(mean( (pred_y - y)^2)) ) %>%
  spread( HOI, RMSE) %>% 
  mutate( RMSE_change = (`1` - `0`)  ) %>% 
  ungroup() %>% 
  mutate( species_lab = factor( species, labels = c('Early', 'Mid', 'Late'))) %>% 
  ggplot( aes( x = species_lab, y = RMSE_change, fill = mod_type)) + 
  geom_bar( stat = 'identity', position = 'dodge') +
  ylab( 'Increase in RMSE') + 
  xlab( 'Species') + 
  my_theme + 
  journal_theme + 
  theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) + 
  guides( color = F)  + 
  ggtitle("A)") + 
  theme(plot.title = element_text(hjust = 0))

RMSE_plot_mod2 <- 
  RMSE_plot_both_models$data %>% 
  filter( mod_type == 'm2') %>%  
  ggplot( aes( x = species_lab, y = RMSE_change, fill = species_lab)) + 
  geom_bar( stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = my_colors[1:3]) + 
  ylab( 'Increase in RMSE') + 
  xlab( 'Species') + 
  guides(fill = F) + 
  journal_theme + 
  ggtitle("A)") + 
  theme(plot.title = element_text(hjust = 0))

error_y_lab <- formula( Average~HOI~effect~(y[obs] - y[pred]))

mean_error_plot_both_mods <- 
  two_sp_df %>% 
  select( -time , -c(X2:X4), -lambda) %>%
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>% 
  filter( n_comp > 0 ) %>% 
  group_by(species, HOI, mod_type) %>% 
  summarise( ME = mean( (y - pred_y) ) ) %>%
  spread( HOI, ME) %>% 
  mutate( ME_change = (`1` - `0`)  ) %>% 
  ungroup() %>% 
  mutate( species_lab = factor( species, labels = c('Early', 'Mid', 'Late'))) %>% 
  mutate( Model = factor(mod_type, labels = c('1', '2'))) %>% 
  ggplot( aes( x = species_lab, y = ME_change, fill = Model)) + 
  geom_bar( stat = 'identity', position = position_dodge()) +
  geom_hline(aes(yintercept = 0)) + 
  scale_color_manual(values = my_colors[1:3])  + 
  ylab( error_y_lab) + 
  xlab( 'Species') + 
  my_theme + 
  journal_theme + 
  theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) + 
  guides( color = F) + 
  ggtitle("B)") + 
  theme(plot.title = element_text(hjust = 0))

mean_error_plot_mod2 <- 
  mean_error_plot_both_mods$data %>% 
  filter( mod_type == 'm2') %>%  
  ggplot( aes( x = species_lab, y = ME_change, fill = species_lab)) + 
  geom_bar( stat = 'identity', position = position_dodge()) +
  geom_hline(aes(yintercept = 0)) + 
  scale_fill_manual(values = my_colors[1:3]) + 
  ylab( error_y_lab) + 
  xlab( 'Species') + 
  journal_theme + 
  guides( fill = F) + 
  ggtitle("B)") + 
  theme(plot.title = element_text(hjust = 0))

error_plots <- grid.arrange(RMSE_plot_mod2, mean_error_plot_mod2, nrow = 1, widths = c(0.49, 0.51))

error_plots_both_mods <- grid.arrange(RMSE_plot_both_models + 
                                        guides( fill = F), 
                                      mean_error_plot_both_mods + 
                                        scale_fill_discrete( 'Model') + 
                                        theme(legend.position = c(0.05, 1), 
                                              legend.background = element_blank(),
                                              legend.justification = c(0,1)), nrow = 1, widths = c(0.49, 0.51))

ggsave( error_plots, filename = 'figures/error_plots.png', width = 7, height = 4)

ggsave( error_plots_both_mods, filename = 'figures/appendix_compare_errors.png', width = 7, height = 4)

