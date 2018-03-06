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
alphas <- matrix( c(1, 0.5, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 20), 3, 3, byrow = T)
lambdas <- c(24, 32, 41)
taus <- c(-1, -1, -0.2)
pars <- list( lambdas = lambdas, alphas = alphas , taus = taus) 

# -------- simulate gradient  -------------------------- # 
maxdens <- 10
base <- 2 
experiments <- expand.grid(N1 = c(0, c(base^c(0:maxdens))), N2 = c(0, base^c(0:maxdens)), N3 = c(0, base^c(0:maxdens)))

# 
form <- as.formula('~ -1 + N1 + N2 + N3')

out <- experiments
for( i in 1:nrow(experiments)){ 
  seeds <- experiments[i, ]
  out[i, ] <- ann_plant_mod(seeds, form, unlist(pars))
}

results <- data.frame(experiments, out)
names(results) <- c(paste0('N', 1:3), paste0('F', 1:3))

results <- 
  results %>%   
  tibble::rownames_to_column('id') %>% 
  filter( (N1 == 0 & N2 == 0) | (N3 == 0 & N2 == 0 ) | (N1 == 0 & N3 == 0 ) ) %>% 
  mutate( lambda =  ifelse(N1 == 0 & N2 == 0 & N3 == 0 , T, F)) %>% 
  gather( competitor, density, N1:N3) %>% 
  filter( lambda | density > 0 ) %>% 
  gather( focal, fecundity, F1:F3)  


results$competitor_label <- paste0( 'competitor\n', results$competitor) 
results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))

fits <- lapply(1:3, function(x, ...) fit_ann_plant(focal = x, ... ), data = results, model = mod_bh, method = 'BFGS' )
converged <- lapply( fits, function(x) x$convergence) == 0
fitpars <- lapply( fits, function(x) x$par )
fits[!converged] <- mapply(x = which(!converged), y = fitpars[!converged], FUN = function(x, y, ... ) fit_ann_plant(focal = x, my_inits = y, ... ), MoreArgs = list( data = results, model = mod_bh, method = 'BFGS' ), SIMPLIFY = F)


fits.2 <- lapply(1:3, function(x, ...) fit_ann_plant(focal = x, ... ), data = results, model = mod_bh2, method = 'L-BFGS-B',  my_inits = c(20, 0,0,0), lower = c(1, 1e-6, 1e-6, 1e-6 ))
converged <- lapply( fits.2, function(x) x$convergence) == 0
fitpars <- lapply( fits.2, function(x) x$par )
fits.2[!converged] <- mapply(x = which(!converged), y = fitpars[!converged], FUN = function(x, y, ... ) fit_ann_plant(focal = x, my_inits = y, ... ), MoreArgs = list( data = results, model = mod_bh2, method = 'BFGS' ), SIMPLIFY = F)

pred_fit1 <- predict_fit(results, model = mod_bh, fits[[1]]$par, 1) 
pred_fit2 <- predict_fit(results, model = mod_bh, fits[[2]]$par, 2)
pred_fit3 <- predict_fit(results, model = mod_bh, fits[[3]]$par, 3)

pred_fit1.2 <- predict_fit(results, model = mod_bh2, fits.2[[1]]$par, 1) 
pred_fit2.2 <- predict_fit(results, model = mod_bh2, fits.2[[2]]$par, 2)
pred_fit3.2 <- predict_fit(results, model = mod_bh2, fits.2[[3]]$par, 3)

figdat <- 
  results %>% 
  left_join(pred_fit1, by = c('id', 'focal')) %>% 
  left_join(pred_fit2, by = c('id', 'focal')) %>% 
  left_join(pred_fit3, by = c('id', 'focal')) %>%
  left_join(pred_fit1.2, by = c('id', 'focal')) %>% 
  left_join(pred_fit2.2, by = c('id', 'focal')) %>% 
  left_join(pred_fit3.2, by = c('id', 'focal')) %>%
  gather( predicted, pred_fecundity, starts_with('pred')) %>% 
  separate( predicted, c('t1', 'model', 'predicted_sp'), sep = '\\.') %>% 
  select(-t1, -predicted_sp) %>% 
  filter( !is.na(pred_fecundity))

all_fits <- 
  figdat %>%
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes( y = pred_fecundity, color = model), linetype = 1) +
  facet_grid(focal_label ~ competitor_label )

ggsave(all_fits, filename = 'figures/fit_ann_plant_model.png', width = 10, height = 5)
ggsave(all_fits + scale_y_continuous(trans = 'log'), filename = 'figures/fit_ann_plant_model_log.png', width = 10, height = 5)

# compare parameters 


fitted <- 
  data.frame( N1 = fits[[1]]$par, N2 = fits[[2]]$par, N3 = fits[[3]]$par) %>%
  mutate( par = c('lambda', paste0('alpha_', 1:3), 'tau')) %>% 
  gather( species, value, N1:N3) %>% 
  mutate( par = str_replace(par, '_', str_extract(species, '\\d+'))) %>% 
  mutate( type = 'fitted')

original <- 
  data.frame( species = c('N1', 'N2', 'N3'), lambda = pars$lambdas, alpha = pars$alphas, tau = pars$taus ) %>% 
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
