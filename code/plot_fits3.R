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
model_labs <- c('Hassel', 'Hassel + HOI', 'Hassel pw', 'Model 2', 'Model 2 + HOI', 'Model 2 pw')  
  
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
  filter( (density != 0)|(m1 == max(m1) ) ) %>%
  gather( model, y_hat, starts_with('m'))  %>% 
  mutate( `Best Fit` = factor(model, labels = model_labs)) %>%  
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
  filter( model %in% c('m1_pw', 'm2_pw')) %>% 
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

top_pw_preds$`Best Fit` <- factor( top_pw_preds$`Best Fit`, labels = model_labs[c(1,4)])

pp1 <- 
  pp + 
  geom_line(data = top_pw_preds %>% distinct(), 
            aes( y = y_hat, linetype = `Best Fit`))  

ggsave(pp1,
       filename = 'figures/figure_3_new.png',
       width = 7,
       height = 4)

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
  mutate( `Best Fit` = factor(model, labels = model_labs)) %>% 
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

legend_pos <- all_res %>% group_by( species ) %>% summarise( ymax = max(y), ymin = min(y))

y_pos <- legend_pos$ymax - 0.05*(legend_pos$ymax - legend_pos$ymin)
y_pos <- c(y_pos, legend_pos$ymax - 0.1*(legend_pos$ymax - legend_pos$ymin))

annotate_df <- data.frame( Species = species_labs2,
                           y_pos = y_pos, 
                           x_pos = c(8, 8, 8), 
                           labels = c('C1 = Mid', 'C1 = Early', 'C1 = Early', 
                                      'C2 = Late', 'C2 = Late', 'C2 = Mid'))


Model1_fits <- two_sp_plot(all_res, 
                           all_preds %>% filter( model %in% c('m1', 'm1_HOI')), 
                           label_df = annotate_df)

Model2_fits <- two_sp_plot(all_res, 
                           all_preds %>% filter( model %in% c('m2', 'm2_HOI')), 
                           label_df = annotate_df)

ggsave(Model1_fits, 
       filename = 'figures/figure_4_new.png', 
       width = 7.5, 
       height = 4.5)

ggsave(Model2_fits, 
       filename = 'figures/figure_S2_new.png', 
       width = 7.5, 
       height = 4.5)

# Plot model error  ---------------------------------------# 

Models <- model_labs[c(1,2,4,5)]
modlist <- list( fit1, fit1HOI, fit2, fit2HOI)
out <- list()
for( i in 1:length(Models)){ 
  mo <- as.numeric( factor( Models))[i]
  temp_error <- do.call( rbind, lapply( modlist[[mo]], function(x) sqrt(  x$m$deviance() )) ) 
  out[[i]] <- data.frame( error = temp_error, Species = species_labs, Model = Models[i])
}


error_theme <- 
  journal_theme + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_text(size = 14), 
        axis.title = element_text(size = 16))

error_df <- 
  do.call(rbind, out ) %>% 
  separate( Model, c('Model', 'HOI'), sep = ' \\+ ', fill = 'right') %>% 
  mutate( HOI = ifelse( is.na(HOI), 'No HOI', 'HOI')) %>% 
  mutate( Model = factor( Model, levels = c('Hassel', 'Model 2'), ordered = T)) %>% 
  mutate( Species = factor( Species, levels = species_labs, ordered = T)) %>% 
  mutate( HOI = factor( HOI, levels = c('No HOI', 'HOI'), ordered = T)) 

gg_error <- 
  error_df %>% 
  ggplot( aes( x = Model, y = error, fill = HOI) ) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    scale_fill_discrete() + 
    facet_wrap(~Species) + 
    ggtitle('Focal species') + 
    error_theme + 
    guides(fill = guide_legend(title = NULL) )  + 
    theme(axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 20))  + 
    ylab( 'Residual Sum-of-Squares')

gg_error

ggsave( gg_error, 
        filename = 'figures/figure_5_new.png', 
        width = 7, 
        height = 4)
  
save(error_df, file = 'output/model_errors.rda')
