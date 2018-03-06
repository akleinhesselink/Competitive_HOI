library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(parallel)

rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

predict_fit <- function( dat, model, pars, foc = 1 ){ 
  
  dat <- 
    dat %>% 
    filter( focal == paste0('F', foc)) %>% 
    distinct(id, focal, competitor, density) %>% 
    spread( competitor, density , fill = 0) %>% 
    arrange(N1, N2, N3) %>% 
    group_by( N1, N2, N3) %>% 
    arrange( as.numeric(id) ) %>% 
    filter( row_number() == 1 )
  
  dat$y <- NA
  dat$pred <- as.numeric(model(pars = pars, dat, form = form, predict = T))
  
  dat %>% 
    ungroup() %>% 
    select(id, focal, pred) %>% 
    filter( !is.na(pred)) %>% 
    distinct() %>% 
    arrange( as.numeric(id)) %>% 
    mutate( focal_predicted = paste0('pred_', focal)) %>% 
    spread( focal_predicted, pred )
}

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

# -------- simulate multiple years -------------------------- # 
maxdens <- 8
experiments <- expand.grid(N1 = c(0, c(2^c(0:maxdens))), N2 = c(0, 2^c(0:maxdens)), N3 = c(0, 2^c(0:maxdens)))

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


fit1 <- fit_ann_plant(results, focal = 1, model = mod_bh, my_inits = c(24, 1, 0, 0, -1), method = 'L-BFGS-B', lower = c(10, 0, 0, 0, -2), upper = c(1000, 20, 20, 20, -0.0001))
fit1 <- fit_ann_plant(results, focal = 1, model = mod_bh, my_inits = fit1$par, method = 'BFGS', lower = c(10, 0, 0, 0, -2), upper = c(1000, 20, 20, 20, -0.0001), control = list(maxit = 1e5))
fit2 <- fit_ann_plant(results, focal = 2, model = mod_bh, method = 'BFGS')

fit3 <- fit_ann_plant(results, focal = 3, model = mod_bh, my_inits = c(41, 0, 0, 10, -0.5), method = 'L-BFGS-B', lower = c(10, 0, 0, 0, -2), upper = c(1000, 20, 20, 100, -0.0001))
fit3 <- fit_ann_plant(results, focal = 3, model = mod_bh, my_inits = fit1$par, method = 'BFGS')


fit1.2 <-fit_ann_plant(results, focal = 1, model = mod_bh2, my_inits = c(41, 0.1, 0.1, 0.1), method = 'BFGS')
fit2.2 <-fit_ann_plant(results, focal = 2, model = mod_bh2, my_inits = c(41, 0.1, 0.1, 0.1), method = 'BFGS')
fit3.2 <-fit_ann_plant(results, focal = 3, model = mod_bh2, my_inits = c(41, 0.1, 0.1, 0.1), method = 'BFGS')
fit3.2

fit3.2 <-fit_ann_plant(results, focal = 3, model = mod_bh2, my_inits = fit3.2$par, method = 'BFGS')
fit3.2$par

pred_fit1 <- predict_fit(results, model = mod_bh, fit1$par, 1) 
pred_fit2 <- predict_fit(results, model = mod_bh, fit2$par, 2)

my_par <- fit3.2$par
pred_fit3 <- predict_fit(results, model = mod_bh2, my_par, 3)

testtest <- 
  results %>% 
  left_join(pred_fit1, by = c('id', 'focal')) %>% 
  left_join(pred_fit2, by = c('id', 'focal')) %>% 
  left_join(pred_fit3, by = c('id', 'focal')) 

testtest <- testtest %>% 
  gather( predicted, pred_fecundity, starts_with('pred')) %>% 
  select(-predicted) %>% 
  filter( !is.na(pred_fecundity))

testtest %>%
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes( y = pred_fecundity), linetype = 2) + 
  facet_grid(focal_label ~ competitor_label )

testtest %>%
  filter( focal == 'F3') %>% 
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes(y = pred_fecundity), alpha = 0.5, linetype = 2) + 
  ylab('fecundity of N3')  + 
  facet_grid( ~ competitor_label )

testtest %>%
  filter( focal == 'F3') %>% 
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes(y = pred_fecundity), alpha = 0.5, linetype = 2) + 
  ylab('fecundity of N3')  + 
  facet_grid( ~ competitor_label ) + 
  scale_y_continuous(trans = 'log')
