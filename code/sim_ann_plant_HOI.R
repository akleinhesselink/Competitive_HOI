
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
betas <- matrix(c(0.1,  0.2,  0.01, 
                  0.001,  0.1,  -0.01, 
                  0.01,  0.02,  0.01), nspp, nspp, byrow = T)

lambdas <- c(24, 32, 41)
taus <- c(-1.01, -1, -0.9)
#pars <- list( lambdas = lambdas, alphas = cbind(alphas, betas), taus = taus ) 
# 
maxdens <- 20
base <- 1.2

# -----------------------------------------------------
experiments <- make_experiments(maxdens, base, nspp)
#experiments <- make_biculture(experiments)

out <- experiments
mm <- model.matrix(formHOI, experiments)

for( i in 1:nspp) { 
  out[,i] <- mod_bh(pars = c(lambdas[i], taus[i], alphas[i, ], betas[i, ]), y = NA, mm = mm, predict = T)
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

par_ests <- expand.grid( focal = 1:nspp, comp = 1:nspp, lambda = NA, tau = NA, alpha = NA)

for( i in 1:nrow(par_ests)) {
  temp <- fit_2_converge(n_seq = 20, 
                         start_sd = 3, 
                         min_sd = 0.2, 
                         init_tau = -1, 
                         init_alpha = 0, 
                         data = results, 
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
                             data = results, 
                             focal = i, 
                             form = form1, 
                             FUN = fit_single_tau)
  
  best_fit <- temp$res[[which.min( unlist( lapply( temp$res, function(x) x$value)))]]
  
  if(best_fit$convergence == 0 ) { 
    fit_pars[[i]] <- best_fit$par 
  }
}

fit_pars <- do.call(rbind, fit_pars)

# now fit HOIs 
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
                             data = results, 
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
round( beta_est, 3) == betas

fits <- data.frame(focal = 1:3, lambda = lambda_est, tau = tau_est)
fits$alpha <- alpha_est
fits$beta <- beta_est

pred_fit1 <- list()
pred_fit2 <- list()

for( i in 1:nrow(fits)){ 
  
  pred_fit1[[i]] <- predict_fit(pars = unlist(fits[i, c('lambda', 'tau', 'alpha')]), 
              foc = i, 
              dat = results, 
              model = mod_bh_ll, 
              form = form1)
  
  pred_fit2[[i]] <- predict_fit(pars = unlist(fits[i, c('lambda', 'tau', 'alpha', 'beta')]), 
                                foc = i, 
                                dat = results, 
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
  left_join(results, figdat, by = c('id', 'focal'))

p1 <- plot_two_sp(all_fits, focal = 'F1', 'N1', 'N2', 'N3') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p2 <- plot_two_sp(all_fits, focal = 'F1', 'N1', 'N3', 'N2') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p3 <- plot_two_sp(all_fits, focal = 'F1', 'N2', 'N3', 'N1') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3)) + 
  theme(legend.position = c(0.85, 0.7) )

sp1_fits <- grid.arrange(p1, p2, p3, ncol = 3)

p1 <- plot_two_sp(all_fits, focal = 'F2', 'N2', 'N1', 'N3') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p2 <- plot_two_sp(all_fits, focal = 'F2', 'N2', 'N3', 'N1') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p3 <- plot_two_sp(all_fits, focal = 'F2', 'N1', 'N3', 'N2') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3)) + 
  theme(legend.position = c(0.85, 0.7) )

sp2_fits <- grid.arrange(p1, p2, p3, ncol = 3)

p1 <- plot_two_sp(all_fits, focal = 'F3', 'N3', 'N1', 'N2') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p2 <- plot_two_sp(all_fits, focal = 'F3', 'N3', 'N2', 'N1') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3), guide = F) + 
  theme(legend.position = c(0.85, 0.7) )

p3 <- plot_two_sp(all_fits, focal = 'F3', 'N1', 'N2', 'N3') + 
  geom_line(aes( y = pred_fecundity, linetype = form)) + 
  scale_linetype_manual(values = c(2,3)) + 
  theme(legend.position = c(0.85, 0.7) )

sp3_fits <- grid.arrange(p1, p2, p3, ncol = 3)

# compare parameters 
fitted <- 
  data.frame( species = paste0('N', 1:nspp), 
              lambda = lambda_est,
              alpha = alpha_est, 
              betas = beta_est, 
              tau = tau_est ) %>%   
  gather( par, value, lambda:tau) %>%
  mutate( par = str_replace(par, '\\.', str_extract(species, '\\d+'))) %>%
  mutate( type = 'fitted')

original <- 
  data.frame( species = paste0('N', 1:nspp), 
              lambda = lambdas, 
              alpha = alphas, 
              betas = betas, 
              tau = taus ) %>% 
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

