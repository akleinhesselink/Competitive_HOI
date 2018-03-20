
rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# set parameters ------------------------------------- 
nspp <- 3 
alphas <- matrix( c(1, 0.5, 0.1,  
                    0.3, 1, 0.4, 
                    0.1, 0.5, 1), nspp, nspp, byrow = T)
#betas <- alphas 
#betas[] <- rep(0.01, nspp*nspp)
betas <- matrix( c(0.01, 0.001, 0.00, 
                   0.00, 0.01, 0.0002, 
                   0.0, 0.0, 0.01), nspp, nspp, byrow = T)

lambdas <- c(24, 32, 41)
taus <- c(-1, -1, -1)
pars <- list( lambdas = lambdas, alphas = cbind(alphas, betas), taus = taus ) 
# 
maxdens <- 10
base <- 2 

# -----------------------------------------------------
experiments <- make_experiments(maxdens, base, nspp)
experiments <- make_monoculture(experiments)

out <- experiments
for( i in 1:nrow(experiments)){ 
  seeds <- experiments[i,1:nspp]
  out[i,1:nspp] <- ann_plant_mod(seeds, formHOI, unlist(pars))
}

names(out)[1:nspp] <- paste0('F', 1:nspp)

results <- 
  left_join(experiments, out, by = 'id') %>% 
  gather( focal, fecundity, starts_with('F')) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  group_by( id, focal ) %>% 
  mutate( comp_n = sum(density > 0)) %>% 
  filter( comp_n == 0 | density > 0 ) %>% 
  spread( competitor, density, fill = 0) %>% 
  ungroup()

results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))


my_fits1 <- fit_both_mods(focal = 1, 
                          n_seq = 20, 
                          start_sd = 2, 
                          min_sd = 0.01, 
                          form1 = form1,
                          data = results, 
                          model = mod_bh_ll, 
                          lower1 = c(1, 0, 0, 0, -2),
                          inits1 = c(10, 0, 0, 0, -1.5), 
                          method = 'L-BFGS-B')

my_fits2 <- fit_both_mods(focal = 2, 
                          n_seq = 20, 
                          start_sd = 2, 
                          min_sd = 0.01, 
                          form1 = form1,
                          data = results, 
                          model = mod_bh_ll, 
                          lower1 = c(10, 0, 0, 0, -2),
                          inits1 = c(10, 0, 0, 0, -1.1), 
                          method = 'L-BFGS-B')

my_fits3 <- fit_both_mods(focal = 3, 
                          n_seq = 20, 
                          start_sd = 2, 
                          min_sd = 0.01, 
                          form1 = form1,
                          data = results, 
                          model = mod_bh_ll, 
                          lower1 = c(10, 0, 0, 0, -2),
                          inits1 = c(10, 0, 0, 0, -1.1), 
                          method = 'L-BFGS-B')

fits1 <- lapply(list(my_fits1, my_fits2, my_fits3), function(x) x$fit1)

pred_fit1 <- mapply( x = 1:nspp, 
                     y = fits1, 
                     z = c(mod_bh_ll, mod_bh_ll, mod_bh_ll), 
                     FUN = function(x,y,z, ... ) 
                       predict_fit(pars = y$par, 
                                   foc = x, 
                                   model = z, 
                                   dat = results, 
                                   form = form1), 
                     SIMPLIFY = F)

preds <- do.call( rbind, lapply( pred_fit1, function(x) x %>% gather( predicted, pred_fecundity, starts_with('pred'))))

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
  geom_line(aes( y = pred_fecundity, color = form), linetype = 1) +
  facet_grid(focal_label ~ competitor + form)

sp3_fits <- 
  all_fits$data %>%
  filter( focal == 'F3') %>% 
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes( y = pred_fecundity, color = form), linetype = 1) +
  ylab('fecundity of N3')  + 
  facet_grid( ~ competitor )

sp3_fits_log <- 
  sp3_fits + scale_y_continuous(trans = 'log')

sp3_fits
sp3_fits_log
all_fits
