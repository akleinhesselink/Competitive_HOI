library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)
rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# set parameters ------------------------------------- 
nspp <- 3 
alphas <- matrix( c(1, 0.5, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 20), nspp, nspp, byrow = T)
lambdas <- c(24, 32, 41)
taus <- c(-1, -1, -0.2)
pars <- list( lambdas = lambdas, alphas = alphas , taus = taus ) 
# 
maxdens <- 10
base <- 2 

# -----------------------------------------------------
experiments <- make_experiments(maxdens , base, nspp)
monocultures <- make_monoculture(experiments)

out <- monocultures
for( i in 1:nrow(monocultures)){ 
  seeds <- monocultures[i,1:nspp]
  out[i,1:nspp] <- ann_plant_mod(seeds, form1, unlist(pars))
}

names(out)[1:nspp] <- paste0('F', 1:nspp)

results <- 
  left_join(monocultures, out, by = 'id') %>% 
  gather( focal, fecundity, starts_with('F')) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  group_by( id, focal ) %>% 
  mutate( comp_n = sum(density > 0)) %>% 
  filter( comp_n == 0 | density > 0 ) %>% 
  spread( competitor, density, fill = 0) %>% 
  ungroup()


results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))

fits.1 <- fit_2_converge(results, model = mod_bh, method = 'L-BFGS-B', form = form1, lower = c(0, 0, 0, 0, -2))
fits.2 <- fit_2_converge(results, model = mod_bh2, method = 'L-BFGS-B', form = form1, my_inits = c(10, 1,1,1), lower = c(1, 0, 0, 0))

pred_fit1 <- mapply( x = 1:nspp, y = fits.1, FUN = function(x,y, ... ) predict_fit(pars = y$par, foc = x, dat = results, model = mod_bh, form = form1), SIMPLIFY = F)
pred_fit2 <- mapply( x = 1:nspp, y = fits.2, FUN = function(x,y, ... ) predict_fit(pars = y$par, foc = x, dat = results, model = mod_bh2, form = form1), SIMPLIFY = F)

preds <- do.call( rbind, lapply( c(pred_fit1, pred_fit2), function(x) x %>% gather( predicted, pred_fecundity, starts_with('pred'))))

figdat <-
  preds %>% 
  separate( predicted, c('t1', 'model', 'predicted_sp'), sep = '\\.') %>% 
  select(-t1, -predicted_sp) %>% 
  filter( !is.na(pred_fecundity))


all_fits <- 
  left_join(results, figdat, by = c('id', 'focal')) %>%
  gather( competitor, density, starts_with('N')) %>% 
  filter( comp_n == 0 | density > 0 ) %>%  
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes( y = pred_fecundity, color = model), linetype = 1) +
  facet_grid(focal_label ~ competitor )

ggsave(all_fits, filename = 'figures/fit_ann_plant_model.png', width = 10, height = 5)
ggsave(all_fits + scale_y_continuous(trans = 'log'), filename = 'figures/fit_ann_plant_model_log.png', width = 10, height = 5)

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

ggsave(pars_plot, filename = 'figures/plot_ann_plant_pars.png')
