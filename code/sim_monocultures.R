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
save(parms, file = 'data/parms.rda')
rm( list = names(parms) )

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# -----------------------------------------------------
nspp <- length(parms$K)
maxdens <- 8
base <- 2 

experiments <- make_experiments(maxdens, base, nspp)

monocultures <- make_monoculture(experiments)

experiments <- 
  rbind( add_focal(monocultures, c(1,0,0), focal_lab = 'F1'), 
       add_focal(monocultures, c(0,1,0), focal_lab = 'F2'), 
       add_focal(monocultures, c(0,0,1), focal_lab = 'F3')) 

results <- experiments
for( i in 1:nrow(results)){ 
  results[i, grep('N', names(results))] <- run_experiment(results[i, grep('N', names(results))], parms)
}

results[ , grep('N', names(results)) ] <- (results %>% select(starts_with('N')))/(experiments %>% select(starts_with('N')))

names( results ) <- str_replace( names(results), '^N', 'F')

results2 <- merge( monocultures, results, by = c('id')) %>% 
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


par_ests <- expand.grid( focal = 1:nspp, comp = 1:nspp, lambda = NA, tau = NA, alpha = NA)

for( i in 1:nrow(par_ests)) {
  temp <- fit_2_converge(n_seq = 20, 
                         start_sd = 3, 
                         min_sd = 0.2, 
                         init_tau = -1, 
                         init_alpha = 0, 
                         data = results2, 
                         focal = par_ests$focal[i], 
                         comp  = par_ests$comp[i], 
                         form = form1)
  
  best_fit <- temp$res[[which.min( unlist( lapply( temp$res, function(x) x$value)))]]
  
  if(best_fit$convergence == 0 ) { 
    par_ests[i, c('lambda', 'tau', 'alpha')] <- best_fit$par 
  }
}

par_ests <- 
  par_ests %>% 
  group_by( focal ) %>% 
  mutate( lambda_hat = mean(lambda), tau_hat = mean(tau)) %>% 
  gather( alpha_par, alpha_value, alpha) %>% 
  unite( alpha, alpha_par, focal, comp, sep = '', remove = F) %>% 
  select( - alpha_par ) %>% 
  arrange( alpha ) 

alpha_est <- matrix( par_ests$alpha_value, nspp, nspp , byrow = T)
lambda_est <-  unique(par_ests$lambda_hat)
tau_est <- unique( par_ests$tau_hat)

fit_pars <- list(NA)
for( i in 1:nspp) {
  temp <- fit_2_converge(n_seq = 20, 
                         start_sd = 3, 
                         min_sd = 0.2, 
                         init_lambda = lambda_est[i],
                         init_tau = tau_est[i], 
                         init_alpha = alpha_est[i, ], 
                         data = results2, 
                         focal = i, 
                         form = form1, 
                         FUN = fit_single_tau)
  
  best_fit <- temp$res[[which.min( unlist( lapply( temp$res, function(x) x$value)))]]
  
  if(best_fit$convergence == 0 ) { 
    fit_pars[[i]] <- best_fit$par 
  }
}
fit_pars <- do.call(rbind, fit_pars)

# now fit HOIs ---------------------------------------------- # 

par_ests <- expand.grid( focal = 1:nspp, comp = 1:nspp)
par_ests$lambda = lambda_est
par_ests$tau = fit_pars[,1]

alpha_est <- fit_pars[, 1 + 1:nspp]
par_ests$alpha <- matrix( rep( t(alpha_est), nspp), nspp*nspp, nspp, byrow = T)
par_ests$beta  <- NA

for(i in 1:nrow(par_ests)){ 
  
  temp <- fit_2_converge(n_seq = 20, 
                         start_sd = 3, 
                         min_sd = 0.2, 
                         init_lambda = par_ests$lambda[i],
                         init_tau = par_ests$tau[i],
                         init_alpha = par_ests$alpha[i,], 
                         data = results2, 
                         focal = par_ests$focal[i], 
                         comp  = par_ests$comp[i], 
                         form = formHOI, 
                         FUN = fit_HOI)
  
  best_fit <- temp$res[[which.min( unlist( lapply( temp$res, function(x) x$value)))]]
  
  if(best_fit$convergence == 0 ) { 
    par_ests$beta[i] <- best_fit$par 
  }
  
}

beta_est <- matrix(par_ests$beta, nspp, nspp)

fits <- data.frame(focal = 1:3, lambda = lambda_est, tau = tau_est)
fits$alpha <- alpha_est
fits$beta <- beta_est

pred_fit1 <- list()
pred_fit2 <- list()

for( i in 1:nrow(fits)){ 
  
  pred_fit1[[i]] <- predict_fit(pars = unlist(fits[i, c('lambda', 'tau', 'alpha')]), 
                                foc = i, 
                                dat = results2, 
                                model = mod_bh_ll, 
                                form = form1)
  
  pred_fit2[[i]] <- predict_fit(pars = unlist(fits[i, c('lambda', 'tau', 'alpha', 'beta')]), 
                                foc = i, 
                                dat = results2, 
                                model = mod_bh_ll, 
                                form = formHOI)
  
  
}

preds <- do.call( rbind, lapply( c(pred_fit1, pred_fit2), function(x) x %>% gather( predicted, pred_fecundity, starts_with('pred'))))

figdat <-
  preds %>% 
  separate( predicted, c('t1', 'model', 'predicted_sp'), sep = '\\.') %>% 
  select(-t1, -predicted_sp) %>% 
  filter( !is.na(pred_fecundity))

all_fits <- 
  left_join(results2, figdat, by = c('id', 'focal')) %>%
  gather( competitor, density, starts_with('N')) %>% 
  filter( comp_n == 0 | density > 0 ) %>%  
  ggplot(aes( x = density, y = fecundity) ) + 
  geom_point() + 
  geom_line(aes( y = pred_fecundity, color = form), linetype = 1) +
  facet_grid(focal_label ~ competitor + form)

all_fits

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

