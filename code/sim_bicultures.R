library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(parallel)
library(stringr)

rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

load('data/parms.rda')

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# -----------------------------------------------------
nspp <- length(parms$K)
maxdens <- 10
base <- 2 

experiments <- make_experiments()

biculture <- make_biculture(experiments)

experiments <- 
  rbind( add_focal(biculture, c(1,0,0), focal_lab = 'F1'), 
         add_focal(biculture, c(0,1,0), focal_lab = 'F2'), 
         add_focal(biculture, c(0,0,1), focal_lab = 'F3')) 

results <- experiments
for( i in 1:nrow(results)){ 
  results[i, grep('N', names(results))] <- run_experiment(results[i, grep('N', names(results))], parms)
}

results[ , grep('N', names(results)) ] <- (results %>% select(starts_with('N')))/(experiments %>% select(starts_with('N')))

names( results ) <- str_replace( names(results), '^N', 'F')

results2 <- merge( biculture, results, by = c('id')) %>% 
  gather( focal2, fecundity, starts_with('F', ignore.case = F)) %>% 
  filter( focal == focal2) %>% 
  select(-focal2) %>%
  gather( competitor, density, starts_with('N')) %>% 
  group_by( id, focal ) %>% 
  mutate( comp_n = sum(density > 0)) %>% 
  spread( competitor, density, fill = 0) %>% 
  ungroup()

results2$focal_label <- paste0( 'focal\n', str_replace( results2$focal, 'F', 'N'))

results2 %>%
  filter( comp_n  < 2 ) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  filter( comp_n == 0 | density > 0 ) %>%
  ggplot(aes( x = density, y = fecundity, color = competitor) ) + 
  geom_point(alpha = 1) + 
  geom_line(alpha = 0.5) + 
  scale_color_discrete(guide = F) + 
  facet_grid(focal_label ~ competitor)

plot_two_sp(results2, focal = 'F1', 'N1', 'N2', 'N3')
plot_two_sp(results2, focal = 'F1', 'N1', 'N3', 'N2')

plot_two_sp(results2, focal = 'F2', 'N2', 'N1', 'N3')
plot_two_sp(results2, focal = 'F2', 'N2', 'N3', 'N1')

plot_two_sp(results2, focal = 'F3', 'N3', 'N1', 'N2')
plot_two_sp(results2, focal = 'F3', 'N3', 'N2', 'N1')

results2 <- 
  results2 %>% 
  filter( N1 + N2 + N3 <= max(N1) + 10 )

fits.1 <- fit_2_converge(results2, model = mod_bh, method = 'BFGS', form = form1)
fits.2 <- fit_2_converge(results2, model = mod_bh2, method = 'L-BFGS-B', form = form1, my_inits = c(10, 1, 1, 1 ), lower = c(1, 1e-30, 1e-30, 1e-30))

fits.3 <- fit_2_converge(results2, model = mod_bh, method = 'BFGS', form = formHOI)



pred_fit1 <- mapply( x = 1:nspp, y = fits.1, FUN = function(x,y, ... ) predict_fit(pars = y$par, foc = x, dat = results2, model = mod_bh, form = form1), SIMPLIFY = F)
pred_fit2 <- mapply( x = 1:nspp, y = fits.2, FUN = function(x,y, ... ) predict_fit(pars = y$par, foc = x, dat = results2, model = mod_bh2, form = form1), SIMPLIFY = F)

preds <- do.call( rbind, lapply( c(pred_fit1, pred_fit2), function(x) x %>% gather( predicted, pred_fecundity, starts_with('pred'))))

figdat <-
  preds %>% 
  separate( predicted, c('t1', 'model', 'predicted_sp'), sep = '\\.') %>% 
  select(-t1, -predicted_sp) %>% 
  filter( !is.na(pred_fecundity))

all_fits <- 
  left_join(results2, figdat, by = c('id', 'focal'))


library(gridExtra)
library(scales)


p1 <- plot_two_sp(all_fits, focal = 'F1', 'N1', 'N2', 'N3') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p2 <- plot_two_sp(all_fits, focal = 'F1', 'N1', 'N3', 'N2') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p3 <- plot_two_sp(all_fits, focal = 'F1', 'N2', 'N3', 'N1') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3)) + 
  theme(legend.position = c(0.85, 0.7) )

sp1_fits <- grid.arrange(p1, p2, p3, ncol = 3)


p1 <- plot_two_sp(all_fits, focal = 'F2', 'N2', 'N1', 'N3') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p2 <- plot_two_sp(all_fits, focal = 'F2', 'N2', 'N3', 'N1') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p3 <- plot_two_sp(all_fits, focal = 'F2', 'N1', 'N3', 'N2') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3)) + 
  theme(legend.position = c(0.85, 0.7) )

sp2_fits <- grid.arrange(p1, p2, p3, ncol = 3)


p1 <- plot_two_sp(all_fits, focal = 'F3', 'N3', 'N1', 'N2') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p2 <- plot_two_sp(all_fits, focal = 'F3', 'N3', 'N2', 'N1') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p3 <- plot_two_sp(all_fits, focal = 'F3', 'N1', 'N2', 'N3') + 
  geom_line(aes( y = pred_fecundity, linetype = model)) + 
  scale_linetype_manual(values = c(2,3)) + 
  theme(legend.position = c(0.85, 0.7) )

sp2_fits <- grid.arrange(p1, p2, p3, ncol = 3)
