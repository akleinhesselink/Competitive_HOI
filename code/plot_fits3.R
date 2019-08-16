rm(list = ls())
source('code/plotting_parameters.R')
load('output/model_fits3.rda')
load('output/processed_results3.rda')

library(tidyverse)
library(gridExtra)
library(grid)

theme1 <- 
  journal_theme + 
  theme( strip.text = element_text(hjust = 0.1), 
         legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5), 
         legend.key.width = unit(2, 'line'),
         legend.text = element_text(size = 14),
         plot.title = element_text(hjust = 0.5, size = 16))

species_labs2 <- paste0( LETTERS[1:3], ') ' , species_labs)

predicted <- 
  predicted %>% 
  ungroup() %>% 
  mutate( id = row_number())  

sim_results <- 
  sim_results %>% 
  ungroup() %>% 
  mutate( id = row_number()) %>% 
  group_by( species ) %>% 
  mutate( lambda_plot = y == max(y)) %>% 
  ungroup()


pw_pred <- 
  predicted %>%   
  filter( n_comp < 2) %>% 
  select( B1:B3, id) %>% 
  gather( comp, density, B1:B3) %>% 
  left_join(predicted, by = 'id') %>%
  filter( (density != 0)|(m1 == max(m1) ) )  %>% 
  gather( model, y_hat, starts_with('m')) %>%
  mutate( `Best Fit` = factor(model)) %>% 
  mutate( Species = factor(species, labels = species_labs2)) %>% 
  mutate( Competitor = factor(comp, labels = species_labs))

pw_results <-
  sim_results %>% 
  filter( n_comp < 2 ) %>% 
  select( B1:B3, id) %>% 
  gather( comp, density, B1:B3) %>%
  filter( density < 15 ) %>% 
  left_join( sim_results, by = 'id') %>% 
  filter( (density != 0)|( lambda_plot ) ) %>% 
  mutate( Species = factor(species, labels = species_labs2)) %>% 
  mutate( Competitor = factor(comp, labels = species_labs)) %>% 
  select( id, comp, density, species, y, n_comp, Species, Competitor, lambda_plot)

top_pw_preds <- 
  pw_pred %>% 
  filter( model %in% c('m1', 'm2')) %>% 
  distinct(Species, density, Competitor, `Best Fit`, y_hat)

pp <- 
  pw_results %>% 
  filter( density > 0 ) %>% 
  ggplot( aes( x = density, y = y, color = Competitor, shape = Competitor)) + 
  geom_point(size = 3) + 
  geom_point(data= pw_results %>% 
               filter( density == 0) %>% distinct(Species, y, density), 
             color = 'black', shape = 1, size = 2, show.legend = F) + 
  facet_wrap(~ Species) + 
  scale_color_manual(values = my_colors[1:3], 
                     name = 'Competitor\nSpecies') +
  scale_shape_manual(values = c(19,17,15), 
                     name = 'Competitor\nSpecies') +  
  scale_linetype_discrete(name = 'Best Fit') + 
  scale_x_continuous(breaks = c(0,4,8)) +
  ylab( 'Per Capita Seed Production') + 
  xlab( 'Competitor Density') + 
  ggtitle('Focal Species') +
  guides(color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
         shape = guide_legend(order = 1)) + 
  theme1



# plot model predictions and "observed" simulated data ---------------- # 
pp1 <- 
  pp + 
  geom_line(data =top_pw_preds, 
            aes( y = y_hat, linetype = `Best Fit`))  
pp1

# ggsave(pp1, 
#        filename = 'figures/figure_3_new.png', 
#        width = 7, 
#        height = 4)

pp2 <- 
  pp + 
  geom_line(data = pw_pred, 
            aes( y = y_hat, linetype = `Best Fit`))  
  
pp2

# ggsave(pw_comp_pred_gg, 
#        filename = 'figures/figure_S1.png', 
#        width = 7, 
#        height = 4)

# Compare predicted effects of two species competition to the observed 
# effects of two species competition. Add separate species effects together 
# to predict total effect. Errors are deviations from additivity of 
# competition. 
theme2 <- 
  journal_theme + 
  theme(strip.text = element_text(hjust = 0.1), 
        legend.background = element_blank())


two_sp_pred <- 
  predicted %>%   
  filter( n_comp < 3, B1 < 15, B2 < 15, B3 < 15) %>% 
  gather( model, y_hat, starts_with('m')) %>%
  mutate( `Best Fit` = factor(model, labels = c('Hassel', 'Hassel + HOI', 'Model 2', 'Model 2 + HOI'))) %>% 
  mutate( Species = factor(species, labels = species_labs2))

two_sp_results <-
  sim_results %>% 
  filter( n_comp < 3, B1 < 15, B2 < 15, B3 < 15 ) %>% 
  mutate( Species = factor(species, labels = species_labs2)) %>% 
  select( id, B1:B3, species, Species, y, n_comp, lambda_plot)

ylab <- 'Per Capita Seed Production'
xlabels <- c('Mid Competitor\nDensity', 'Early Competitor\nDensity', 'Early Competitor\nDensity')
textlabel <- c('Late Competitor\nDensity','Late Competitor\nDensity', 'Mid Competitor\nDensity')

make_2comp_df <- function(x, i ) { 
  cur_sp <- paste0('Y', i)
  intra_comp <- paste0('B', i)
  inter_comp <- paste0('B', c(1:3)[-i])

  x %>% 
    filter( species == cur_sp, get(intra_comp) == 0  ) %>% 
    select(   - !!(intra_comp)) %>%  
    gather( comp_label, density, inter_comp)  %>% 
    mutate( comp_label = 
              factor(comp_label, 
                     labels = c('Competitor 1', 
                                'Competitor 2'))) %>% 
    spread( comp_label, density ) %>% 
    filter( `Competitor 1` < 15) %>% 
    filter( `Competitor 2` %in% c(0, 2, 8)) %>%
    mutate( `Competitor 2` = factor(`Competitor 2`))    
}


all_preds <- 
  do.call( bind_rows, lapply(1:3, FUN = function(i) make_2comp_df(x = two_sp_pred, i )) )

all_res <- 
  do.call( bind_rows, lapply(1:3, FUN = function(i) make_2comp_df(x = two_sp_results, i)))


two_sp_plot <- function( d1, d2, label_df ){ 
  d1 %>% 
    ggplot(aes( x = `Competitor 1`, y = y)) + 
    geom_point(aes( color = `Competitor 2`, shape = `Competitor 2`), size= 3) + 
    geom_line(data = d2, aes( y = y_hat, linetype = `Best Fit`, color = `Competitor 2` )) + 
    facet_wrap(~Species, scales = 'free') + 
    scale_color_manual(values = c('black', 'red', 'blue')) + 
    scale_shape_manual(values = c(1, 0, 2)) + 
    scale_x_continuous(breaks = c(0,4,8), limits = c(0, 8)) + 
    ylab('Per Capita Seed Production') + 
    xlab('Density of Competitor 1') + 
    ggtitle('Focal Species') + 
    guides( 
      color = guide_legend(order = 1, 
                           override.aes = list(linetype = 0), 
                           title = 'Density of\nCompetitor 2'), 
      shape = guide_legend(order = 1, title = 'Density of\nCompetitor 2'), 
      linetype = guide_legend(title = 'Fit')) + 
    theme2 + 
    theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
    geom_text(data = label_df, 
              aes( x = x_pos, y = y_pos, label = labels), size = 4, hjust = 1, show.legend = F)
    
}

annotate_df <- data.frame( Species = species_labs2,
                           y_pos = c(14, 20, 33, 13.9, 19.1, 31.4), 
                           x_pos = c(8, 8, 8), 
                           labels = c('C1 = Mid', 'C1 = Early', 'C1 = Early', 
                                      'C2 = Late', 'C2 = Late', 'C2 = Mid'))


Model1_fits <- two_sp_plot(all_res, 
                           all_preds %>% filter( model %in% c('m1', 'm1_HOI')), 
                           label_df = annotate_df)

Model2_fits <- two_sp_plot(all_res, 
                           all_preds %>% filter( model %in% c('m2', 'm2_HOI')), 
                           label_df = annotate_df)

Model1_fits

ggsave(Model1_fits, 
       filename = 'figures/figure_4_new.png', 
       width = 7.5, 
       height = 4.5)


# Plot average error in two species communities ---------------------------------------# 
error <- 
  sim_results %>% 
  select( B1:B3, species, y) %>% 
  left_join(
    predicted %>% 
      gather( model, y_hat, starts_with('m')) %>% 
      select( model, y_hat, B1:B3, species, n_comp)) %>% 
  filter( !is.na(model)) %>% 
  mutate( deviation = y - y_hat, SE = deviation^2 ) %>% 
  group_by( species, model) %>% 
  summarise ( MSE = mean( SE ),  ME = mean(deviation))  


error_data <- 
  sim_results %>% 
  select( B1:B3 , species, y )  %>% 
  left_join(
    predicted %>% 
      gather( model, y_hat, starts_with('m')) %>% 
      select( model, y_hat, B1:B3, species, n_comp)) %>% 
  filter( !is.na(model)) %>% 
  mutate( HOI = as.numeric( str_detect(model, 'HOI') )) %>% 
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>% 
  filter( n_comp > 0 ) %>% 
  group_by(species, model, HOI) %>%
  summarise( RMSE = sqrt(mean( (y - y_hat)^2)), ME = mean (y - y_hat) ) %>%
  gather( errortype, error, RMSE, ME)   %>% 
  spread( HOI, error)  
  mutate( error_change = (`1` - `0`)  ) %>% 
  ungroup()  %>% 
  mutate( species_lab = factor( species, labels = c('Early', 'Mid', 'Late'))) %>% 
  mutate( Model = factor(mod_type, labels = c('1', '2')))

  
error_theme <- 
  journal_theme + 
  theme(plot.title = element_text(hjust = 0), 
        axis.text.x = element_text(size = 14))

error_y_lab <- formula( Average~HOI~effect~(y[obs] - y[pred]))

RMSE_plot_both_models <- 
  error_data %>% 
  filter( errortype == 'RMSE' ) %>% 
  ggplot( aes( x = species_lab, y = error_change, fill = mod_type)) + 
  geom_bar( stat = 'identity', position = 'dodge') +
  geom_hline(aes(yintercept = 0)) + 
  ylab( 'Increase in RMSE') + 
  xlab( 'Species') + 
  ggtitle("A)") + 
  error_theme

RMSE_plot_mod2 <- 
  error_data %>% 
  filter( errortype == 'RMSE', mod_type == 'm2') %>%  
  ggplot( aes( x = species_lab, y = error_change, fill = species_lab)) + 
  geom_bar( stat = 'identity', position = 'dodge') +
  geom_hline(aes(yintercept = 0)) + 
  scale_fill_manual(values = my_colors[1:3]) + 
  ylab( 'Increase in RMSE') + 
  xlab( 'Species') + 
  guides(fill = F) + 
  ggtitle("A)") + 
  error_theme

ME_plot_both_mods <- 
  error_data %>% 
  filter( errortype == 'ME' ) %>% 
  ggplot( aes( x = species_lab, y = error_change, fill = Model)) + 
  geom_bar( stat = 'identity', position = position_dodge()) +
  geom_hline(aes(yintercept = 0)) + 
  scale_fill_discrete( 'Model') +
  ylab( error_y_lab) + 
  xlab( 'Species') + 
  ggtitle("B)") + 
  error_theme + 
  theme(legend.position = c(0.05, 1), 
        legend.background = element_blank(),
        legend.justification = c(0,1))


ME_plot_mod2 <- 
  error_data %>% 
  filter( mod_type == 'm2', errortype == 'ME' ) %>%  
  ggplot( aes( x = species_lab, y = error_change, fill = species_lab)) + 
  geom_bar( stat = 'identity', position = position_dodge()) +
  geom_hline(aes(yintercept = 0)) + 
  scale_fill_manual(values = my_colors[1:3]) + 
  ylab( error_y_lab) + 
  xlab( 'Species') + 
  guides( fill = F) + 
  ggtitle("B)") + 
  error_theme

error_plots <- 
  grid.arrange(RMSE_plot_mod2, 
               ME_plot_mod2, 
               nrow = 1, 
               widths = c(0.49, 0.51))

error_plots_both_mods <- 
  grid.arrange(RMSE_plot_both_models + 
                 guides( fill = F), 
               ME_plot_both_mods, 
               nrow = 1, 
               widths = c(0.49, 0.51))

ggsave( error_plots, 
        filename = 'figures/figure_5.png', 
        width = 7, 
        height = 4)

ggsave( error_plots_both_mods, 
        filename = 'figures/figure_S2.png', 
        width = 7, 
        height = 4)