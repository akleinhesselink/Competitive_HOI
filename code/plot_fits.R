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
  theme( strip.text = element_text(hjust = 0.5), 
         legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5), 
         legend.key.width = unit(2, 'line'),
         legend.text = element_text(size = 14),
         plot.title = element_text(hjust = 0.5, size = 20), 
         axis.text = element_text(size = 14), 
         axis.title = element_text(size = 16))

xbreaks <- c(0,10,20,30,40)

model_labs <- c('Hassel', 'HOI', 'Model 2')  
  
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
  filter( model == 'm1_pw' )  %>% 
  mutate( `Best Fit` = factor(model, labels = model_labs[1])) %>%  
  mutate( Species = factor(species, labels = species_labs)) %>% 
  mutate( Competitor = factor(comp, labels = species_labs))

pw_results <-
  sim_results %>% 
  filter( n_comp < 2 ) %>% 
  select( B1:B3, id) %>% 
  gather( comp, density, B1:B3) %>%
  left_join( sim_results, by = 'id') %>% 
  filter( (density != 0)|( lambda_plot ) ) %>% 
  mutate( Species = factor(species, labels = species_labs)) %>% 
  mutate( Competitor = factor(comp, labels = species_labs)) %>% 
  select( id, comp, density, species, y, n_comp, Species, Competitor, lambda_plot) 
  
top_pw_preds <- 
  pw_pred %>% 
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
  scale_linetype_discrete(name = 'Model Fit') + 
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(limits = c(10, 160)) + 
  ylab( 'Per Capita Seed Production') + 
  xlab( 'Competitor Density') + 
  ggtitle('Focal Species') +
  guides(color = guide_legend(order = 1, override.aes = list(linetype = 0)), 
         shape = guide_legend(order = 1)) + 
  theme1


# plot model predictions and "observed" simulated data ---------------- # 

top_pw_preds$`Best Fit` <- factor( top_pw_preds$`Best Fit`, labels = model_labs[c(1)])

pp1 <- 
  pp + 
  geom_line(data = top_pw_preds %>% distinct(), 
            aes( y = y_hat, linetype = `Best Fit`)) + 
  geom_text( data = top_pw_preds %>% 
                   distinct(Species) %>% 
                   mutate( Competitor = 'Early',
                           label = paste0(LETTERS[1:3], ')'), 
                           x_pos = -Inf, y_pos = Inf), 
                 aes( x = x_pos, 
                      y = y_pos, 
                      label = label), 
                 show.legend = F, 
                 hjust = -0.5, 
                 vjust = 1.5, 
                 color = 'black', 
                 size = 5)


pp1
ggsave(pp1,
       filename = 'figures/example_fits_fig5.png',
       width = 8,
       height = 5)

# Plot model error  ---------------------------------------# 
model_labs
Models <- model_labs[c(1,2,3)]

modlist <- list( fit1, fit1HOI, fit2)
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
  mutate( Species = factor( Species, levels = species_labs, ordered = T)) 

error_df

gg_error <- 
  error_df %>% 
  ggplot( aes( x = Model, y = error, fill = Species) ) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    facet_wrap(~Species) + 
    ggtitle('Focal species') + 
    guides(fill = "none") + 
    scale_fill_manual(values = my_colors[1:3] ) + 
    error_theme + 
    theme(axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 20))  + 
    ylab( 'Residual Sum-of-Squares') 

gg_error

ggsave( gg_error, 
        filename = 'figures/model_error_barplot.png', 
        width = 8, 
        height = 5)
  
save(error_df, file = 'output/model_errors.rda')

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
  filter( n_comp < 3, B1 < 40, B2 < 40, B3 < 40) %>% 
  gather( model, y_hat, starts_with('m')) %>% 
  filter( model != 'm1_pw') %>% 
  mutate( `Best Fit` = factor(model, labels = model_labs)) %>% 
  mutate( Species = factor(species, labels = species_labs))

two_sp_results <-
  sim_results %>% 
  filter( n_comp < 3, B1 < 40, B2 < 40, B3 < 40 ) %>% 
  mutate( Species = factor(species, labels = species_labs)) %>% 
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
    filter( `Competitor 1` < 40) %>% 
    filter( `Competitor 2` %in% c(0, 1, 16)) %>%
    mutate( `Competitor 2` = factor(`Competitor 2`))    
}


all_preds <- 
  do.call( bind_rows, lapply(1:3, FUN = function(i) make_2comp_df(x = two_sp_pred, i )) )

all_res <- 
  do.call( bind_rows, lapply(1:3, FUN = function(i) make_2comp_df(x = two_sp_results, i)))


two_sp_plot <- function( d1, d2, label_df ){ 
  c1_label <- 'Density of Competitor 1'
  c2_label <- 'Density of\nComp. 2'
  d1 %>% 
    ggplot(aes( x = `Competitor 1`, y = y)) + 
    geom_point(aes( color = `Competitor 2`, shape = `Competitor 2`), size= 3) + 
    geom_line(data = d2, aes( y = y_hat, linetype = `Best Fit`, color = `Competitor 2` )) + 
    facet_wrap(~Species, scales = 'free') + 
    scale_color_manual(values = c('black', 'red', 'blue')) + 
    scale_shape_manual(values = c(1, 0, 2)) + 
    scale_x_continuous(breaks = xbreaks) + 
    ylab('Per Capita Seed Production') + 
    xlab(c1_label) + 
    ggtitle('Focal Species') + 
    guides( 
      color = guide_legend(order = 1, 
                           override.aes = list(linetype = 0), 
                           title = c2_label), 
      shape = guide_legend(order = 1, title = c2_label), 
      linetype = guide_legend(title = 'Best Fit')) + 
    theme2 + 
    theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
    geom_text(data = label_df, 
              aes( x = x_pos, y = y_pos, label = labels), size = 4, hjust = 1, show.legend = F)
  
}

legend_pos <- all_res %>% group_by( species ) %>% summarise( ymax = max(y), ymin = min(y))

y_pos <- legend_pos$ymax - 0.05*(legend_pos$ymax - legend_pos$ymin)
y_pos <- c(y_pos, legend_pos$ymax - 0.12*(legend_pos$ymax - legend_pos$ymin))

annotate_df <- data.frame( Species = species_labs,
                           y_pos = y_pos, 
                           x_pos = c(40, 40, 40), 
                           labels = c('C1=Mid', 'C1=Early', 'C1=Early', 
                                      'C2=Late', 'C2=Late', 'C2=Mid'))


error2 <- error_df %>% 
  mutate( Species = factor(Species, labels = species_labs, ordered = F)) %>% 
  mutate( RSS = round( error, 2)) %>% 
  mutate( my_expression = paste0( 'RSS[', str_replace(Model, ' ', '~'), ']', '==', RSS )) 

error2

y_pos3 <- c(legend_pos$ymax - 0.25*(legend_pos$ymax - legend_pos$ymin))
y_pos4 <- c(legend_pos$ymax - 0.32*(legend_pos$ymax - legend_pos$ymin))

y_posRSS <- c(y_pos3, y_pos4, y_pos3)
error2

error2$x_pos <- 40 
error2$y_pos <- y_posRSS

Model1_fits <- two_sp_plot(all_res, 
                           all_preds %>% filter( model %in% c('m1', 'm1_HOI')), 
                           label_df = annotate_df)


Model2_fits <- two_sp_plot(all_res, 
                           all_preds %>% filter( model %in% c('m2', 'm1_HOI')), 
                           label_df = annotate_df)

Model1_fits <- Model1_fits + 
  geom_text( data = error2 %>% filter( Model != 'Model 2'), 
             aes( x = x_pos, y = y_pos, 
                  label = my_expression), hjust = 1, parse = T)

Model2_fits <- Model2_fits + 
  geom_text( data = error2 %>% filter( Model != 'Hassel'), 
             aes( x = x_pos, y = y_pos, 
                  label = my_expression), hjust = 1, parse = T)


ggsave(Model1_fits, 
       filename = 'figures/model1_fits_supporting_info.png', 
       width = 8, 
       height = 5)

ggsave(Model2_fits, 
       filename = 'figures/model2_fits_supporting_info.png', 
       width = 8, 
       height = 5)


# plot 1:1 
res <- 
  sim_results %>% 
  arrange(species)  %>% 
  select( starts_with('B'), species, y )  %>% 
  split( .$species )


m1 <- mapply(x = fit1, y =res, function(x,y) 1/exp(predict(x, newdata = y)), SIMPLIFY = F) 
m2 <- mapply(x = fit2, y = res, function(x,y) 1/exp(predict(x, newdata = y)), SIMPLIFY = F) 
mHOI <- mapply(x = fit1HOI, y = res, function(x,y) 1/exp(predict(x, newdata = y)), SIMPLIFY = F) 

res <- do.call(bind_rows, res )

res$m1 <- unlist( m1 )
res$m2 <- unlist( m2 )
res$mHOI <- unlist( mHOI)

res <- 
  res %>% 
  gather( model, y_hat, starts_with('m')) 

res$Model <- factor(res$model, labels = model_labs[c(1,3,2)])
res$Species <- factor( res$species, labels = species_labs)

error_df2 <- 
  res %>% 
  mutate( dev = (res$y_hat - res$y)^2 ) %>% 
  group_by( Species, Model ) %>% 
  summarise( RMSE = sqrt( sum(dev)/n() ) )

error_df2 <- 
  error_df2 %>% 
  left_join(
    res %>% 
    ungroup %>% 
    mutate( ymx = max(y)) %>% 
    distinct(Species, ymx) %>% 
    mutate( pos1 = ymx - ymx*0.85, 
            pos2 = ymx - ymx*0.0), 
    by = 'Species') %>% 
  mutate( RMSE = round( RMSE, 2)) %>% 
  mutate( my_expression = paste0( 'italic(RMSE', '==', RMSE, ')' )) 


one2one_plot <- 
  res %>%  
  ggplot( aes( x = y_hat, y = y, color = Species)) + 
  geom_point() + 
  geom_abline(aes( intercept = 0 , slope = 1)) + 
  geom_text( data = error_df2, aes( x = pos2, y = pos1, label = my_expression), hjust = 1, parse = T, size = 4.5, color = 1) +
  facet_grid( Model ~ Species) + 
  xlab('Predicted') + 
  ylab('Observed') + 
  ggtitle('Focal Species') + 
  scale_color_manual(values = my_colors[1:3], guide = F) + 
  theme(plot.title = element_text( hjust = 0.5)) + 
  error_theme

one2one_plot <- 
  one2one_plot +
  geom_text( data = res %>% 
  distinct(Model, Species ) %>% 
  mutate( label = paste0( LETTERS[1:9], ')' ), 
          x_pos = -Inf, 
          y_pos = Inf ), 
  aes( x = x_pos, y = y_pos, label = label), 
  hjust = -0.5, 
  vjust = 1.5 , 
  show.legend = F, color = 1)

one2one_plot %>% 
  ggsave( filename = 'figures/one_to_one_error_plot_fig6.png', 
          width = 8, 
          height = 5)
