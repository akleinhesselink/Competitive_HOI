library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)
library(gridExtra)
library(scales)

rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# set parameters ------------------------------------- 
nspp <- 3 
alphas <- matrix( c(1, 0.5, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), nspp, nspp, byrow = T)
lambdas <- c(24, 32, 41)
taus <- c(-1, -1, -0.5)
pars <- list( lambdas = lambdas, alphas = alphas , taus = taus ) 
# 
maxdens <- 10
base <- 2 

# -----------------------------------------------------
experiments <- make_experiments(maxdens , base, nspp)
biculture <- make_biculture(experiments)

out <- biculture
for( i in 1:nrow(biculture)){ 
  seeds <- biculture[i,1:nspp]
  out[i,1:nspp] <- ann_plant_mod(seeds, form1, unlist(pars))
}

names(out)[1:nspp] <- paste0('F', 1:nspp)

results <- 
  left_join(biculture, out, by = 'id') %>% 
  gather( focal, fecundity, starts_with('F')) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  group_by( id, focal ) %>% 
  mutate( comp_n = sum(density > 0)) %>% 
  filter( comp_n == 0 | density > 0 ) %>% 
  spread( competitor, density, fill = 0) %>% 
  ungroup()


results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))

fits.1 <- fit_2_converge(results, model = mod_bh, method = 'L-BFGS-B', form = form1, lower = c(0, 0, 0, 0, -2))
fits.2 <- fit_2_converge(results, model = mod_bh2, method = 'L-BFGS-B', form = form1, my_inits = c(10, 1,1,1), lower = c(1, 1e-30, 1e-30, 1e-30))

pred_fit1 <- mapply( x = 1:nspp, y = fits.1, FUN = function(x,y, ... ) predict_fit(pars = y$par, foc = x, dat = results, model = mod_bh, form = form1), SIMPLIFY = F)
pred_fit2 <- mapply( x = 1:nspp, y = fits.2, FUN = function(x,y, ... ) predict_fit(pars = y$par, foc = x, dat = results, model = mod_bh2, form = form1), SIMPLIFY = F)

preds <- do.call( rbind, lapply( c(pred_fit1, pred_fit2), function(x) x %>% gather( predicted, pred_fecundity, starts_with('pred'))))

figdat <-
  preds %>% 
  separate( predicted, c('t1', 'model', 'predicted_sp'), sep = '\\.') %>% 
  select(-t1, -predicted_sp) %>% 
  filter( !is.na(pred_fecundity))

all_fits <- 
  left_join(results, figdat, by = c('id', 'focal'))

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

# compare parameters 

fitted <- 
  data.frame( N1 = fits.1[[1]]$par, N2 = fits.1[[2]]$par, N3 = fits.1[[3]]$par) %>%
  mutate( par = c('lambda', paste0('alpha_', 1:nspp), 'tau')) %>% 
  gather( species, value, starts_with('N')) %>% 
  mutate( par = str_replace(par, '_', str_extract(species, '\\d+'))) %>% 
  mutate( type = 'fitted')

original <- 
  data.frame( species = paste0('N', 1:nspp), lambda = pars$lambdas, alpha = pars$alphas, tau = pars$taus ) %>% 
  gather( par, value, lambda:tau) %>%
  mutate( par = str_replace(par, '\\.', str_extract(species, '\\d+'))) %>%
  mutate( type = 'original')

pars_df <- 
  bind_rows(fitted, original ) %>% 
  mutate( par_type = str_extract(par, '[a-z]+')) %>% 
  mutate( par = ifelse(!str_detect(par, '\\d+'), paste0(par, str_extract(species, '\\d+')), par))

pars_plot <- 
  ggplot(pars_df, aes( x = par, y = value, shape = type, color = type)) + 
  geom_point(alpha = 1, size = 4)  +
  scale_shape_manual(values = c(3,1)) + 
  facet_wrap( ~ par_type , scales = 'free') + 
  coord_flip()


pars_plot

