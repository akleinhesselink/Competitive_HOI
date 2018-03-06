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

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.9, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(98, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedlings <- c(1, 0, 0)      # number of seedlings 
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)
rm( list = names(parms) )

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# -------- simulate gradient  -------------------------- # 
maxdens <- 10
base <- 2 
experiments <- expand.grid(N1 = c(0, c(base^c(0:maxdens))), N2 = c(0, base^c(0:maxdens)), N3 = c(0, base^c(0:maxdens)))

form <- as.formula('~ -1 + N1 + N2 + N3')

monocultures <- 
  experiments %>% 
  filter( (N1 == 0 & N2 == 0) | (N3 == 0 & N2 == 0 ) | (N1 == 0 & N3 == 0 ) ) %>% 
  filter( N1 + N2 + N3 > 0)

bicultures <- 
  rbind( 
  experiments %>% 
  filter( N1 == 1 & (N2 == 0 | N3 == 0)) , 
  experiments %>% 
  filter( N2 == 1 & (N1 == 0 | N3 == 0)) , 
  experiments %>% 
  filter( N3 == 1 & (N1 == 0 | N2 == 0) ))
  
experiments <- 
  rbind(monocultures, bicultures) %>% 
  tibble::rownames_to_column('id') %>% 
  gather( competitor, density, N1:N3) %>% 
  group_by( id ) %>%  
  mutate( lambda = ifelse( sum(density) == 1, T, F))  %>% 
  spread( competitor, density) %>% 
  ungroup() 

results <- experiments %>% select(starts_with('N'))

for( i in 1:nrow(experiments)){ 
  results[i, ] <- run_experiment(results[i,], parms)
}

results <- results/(experiments %>% select(starts_with('N')))
names( results ) <- paste0('F', 1:ncol(results))

results <- 
  data.frame(experiments, results ) %>% 
  arrange(as.numeric(id)) %>% 
  gather( focal, fecundity, F1:F3) %>%
  gather( competitor, density, N1:N3) %>%  
  filter( !is.na(fecundity)) %>%
  mutate( density = ifelse( str_extract(competitor, '\\d+') == str_extract(focal, '\\d+'), density - 1, density)) %>% 
  filter( !is.na(density), !is.na(fecundity)) %>% 
  filter( (density > 1) | lambda ) %>% 
  group_by(focal, competitor, density ) %>% 
  arrange(id) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

results$competitor_label <- paste0( 'competitor\n', results$competitor) 
results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))

results %>%
  ggplot(aes( x = density, y = fecundity, color = competitor_label) ) + 
  geom_point(alpha = 1) + 
  geom_line(alpha = 0.5) + 
  facet_grid(focal_label ~ competitor_label )

form <- as.formula('~ -1 + N1 + N2 + N3')

fits <- lapply(1:3, function(x, ...) fit_ann_plant(focal = x, ... ), data = results, model = mod_bh, method = 'BFGS' )
converged <- lapply( fits, function(x) x$convergence) == 0
fitpars <- lapply( fits, function(x) x$par )
fits[!converged] <- mapply(x = which(!converged), y = fitpars[!converged], FUN = function(x, y, ... ) fit_ann_plant(focal = x, my_inits = y, ... ), MoreArgs = list( data = results, model = mod_bh, method = 'BFGS', control = list( maxit = 1e4, reltol = 1e-10)), SIMPLIFY = F)

fit_ann_plant(focal = 3, model = mod_bh, my_inits = c(40, 0, 0, 0, -0.5), data = results, method = 'L-BFGS-B', lower = c(0, 0,0,0,-4), upper= c(1e5, 1e5, 1e5, 1e6, 0))
fits[[3]] <- fit_ann_plant(focal = 3, model = mod_bh, my_inits = c(400, 200, 200, 200, -0.5), data = results, method = 'L-BFGS-B', lower = c(41, 0,0,0, -4 ), upper= c(50, 1e5, 1e5, 1e6, 0))

fits.2 <- lapply(1:3, function(x, ...) fit_ann_plant(focal = x, ... ), data = results, model = mod_bh2, method = 'BFGS',  my_inits = c(40, 1,1,1))
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


sp3_fits <- 
  figdat %>%
  filter( focal == 'F3') %>% 
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes( y = pred_fecundity, color = model), linetype = 1) +
  ylab('fecundity of N3')  + 
  facet_grid( ~ competitor_label )


sp3_fits_log <- 
  sp3_fits + scale_y_continuous(trans = 'log')


ggsave(all_fits, filename = 'figures/mechanistic_model_fits.png', width = 10, height = 5)
ggsave(sp3_fits, filename = 'figures/mechanistic_sp_3_fits.png', width = 10, height = 5)
ggsave(sp3_fits_log, filename = 'figures/mechanisitic_sp_3_fits_log.png', width = 10, height = 5)

