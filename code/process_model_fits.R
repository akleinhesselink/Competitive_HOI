library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)

rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

load('data/rs.rda')
load('data/rs_results.rda')
load('data/model_fits.rda')

rs <- do.call(rbind, rs)
seeds <- do.call(rbind, rs_results)
lambda <- apply( seeds[ apply(rs, 1, sum) == 1 , ], 2, max)

errors <- lapply( model_fits, function(x) lapply( x, function(y) y$fit$value ))

model_errors <- 
  data.frame( error  = unlist(errors, recursive = T) ) %>% 
  tibble::rownames_to_column('id') %>% 
  separate(id, c('species', 'model'), sep = '\\.') 

model_errors

ggplot( model_errors, aes( x = model , y = error)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~ species, scales = 'free')


parms <- unlist( lapply( model_fits, function(x) lapply( x, function(y) y$fit$par)), recursive = F)

l_p <- rep(1, length(parms))
a_p <- c(rep(list(2:4), length(parms)))
b_p <- c(rep(list(NA, c(5:7), NA, c(5:7)), 3))
t_p <- c(rep(list(5, 8, 5:7, 8:13), 3))


fit_lambdas <- mapply( x  = parms, y = l_p, FUN = function(x,y) x[y], SIMPLIFY = T)
fit_alphas <- mapply( x  = parms, y = a_p, FUN = function(x,y) x[y], SIMPLIFY = F)
fit_betas <- mapply( x  = parms, y = b_p, FUN = function(x,y) x[y], SIMPLIFY = F)
fit_taus <- mapply( x  = parms, y = t_p, FUN = function(x,y) x[y], SIMPLIFY = F)


alpha_df <- 
  as.data.frame(t(data.frame(fit_alphas))) %>% 
  tibble::rownames_to_column('id') %>%
  separate( id, c('species','model'), '\\.')  %>% 
  gather( par, val, Var2:Var4) %>% 
  mutate( par = factor( par, labels = c('1', '2', '3'))) %>% 
  mutate( par = paste0( 'a_', str_extract(species, '\\d'), par ))


beta_df <- 
  as.data.frame( t( data.frame( fit_betas[ str_detect( names(fit_betas), 'HOI') ] ))) %>% 
  tibble::rownames_to_column('id') %>%
  separate( id, c('species','model'), '\\.')  %>% 
  gather( par, val, Var5:Var7) %>% 
  mutate( par = factor( par, labels = c('(1x2)', '(1x3)', '(2x3)'))) %>% 
  mutate( par = paste0( 'b_', str_extract(species, '\\d'), par ))

lambda_df <- 
  data.frame( val  = fit_lambdas ) %>% 
  tibble::rownames_to_column('id') %>%
  separate( id, c('species','model'), '\\.') %>% 
  mutate( par = paste0( 'lam_', str_extract(species, '\\d')))

tau_df <- data.frame ( val = unlist( fit_taus[ str_detect( names(fit_taus), 'm2') ] , recursive = F) ) %>% 
  tibble::rownames_to_column('id') %>%
  separate( id, c('species','model', 'par'), '\\.')  

tau_df$par <-  c('1', '2', '3', '1', '2', '3', '(1x2)', '(1x3)', '(2x3)')

tau_df <- 
  tau_df %>% 
  mutate( par = paste0( 'tau_', str_extract(species, '\\d'), par ))

fit_pars_df <- rbind( tau_df, beta_df, alpha_df, lambda_df )


par_plot  <- 
  fit_pars_df %>% 
  mutate( par_type = str_extract(par, '[a-z]+')) %>%
  group_by( par_type) %>% 
  do( gg = ggplot(data = .,  aes( x = par, y = val, fill = species)) + 
    geom_bar(stat = 'identity', position = 'dodge') + facet_wrap( ~ model) )


par_plot$gg[[1]] + coord_flip()
par_plot$gg[[2]] + coord_flip()
par_plot$gg[[3]] + coord_flip()
par_plot$gg[[4]] + coord_flip()



pred_dat <- do.call( rbind, unlist( lapply( model_fits, function(x) lapply( x, function(y) y$data ) ), recursive = F) )

pred_dat <- 
  pred_dat %>% 
  tibble::rownames_to_column('id') %>% 
  separate(id, c('species', 'model', 'row'), sep = '\\.' ) 

pred_dat %>% head

test2 %>% 
  filter( one_species) %>%
  gather( comp, num , N1:N3 ) %>% 
  filter( num_sp == 0 | num > 0) %>% 
  ggplot( aes( x = num, y = y, color = comp)) + 
  geom_point() + 
  geom_line(aes( y = pred), linetype = 2)





#
test <- results
names(test) <-  c('N1', 'N2', 'N3')
test <- lapply( test, function( x) {names(x) <- c('m1', 'm2'); x} )
result_dat <- data.frame( value = unlist( test), label = names(unlist(test)))

result_dat <- 
  result_dat %>% 
  separate( label, c('species', 'model', 'par'), sep = '\\.', extra = 'warn', fill = 'left', convert = T) %>% 
  filter( par != 'counts') %>% 
  group_by(species, model) %>% 
  mutate( type = str_extract( par , '[a-z]+')) %>% 
  mutate( par = ifelse( type == 'par' & par == 'par', 'lambda', par)) %>%
  group_by( species, model, type ) %>% 
  arrange(par) %>% 
  mutate(par = ifelse( type == 'par'  & row_number() == max(row_number()), 'tau', par)) %>% 
  mutate( par = str_replace(par, 'par', 'alpha')) %>% 
  arrange( species, model, par)

save(results, file = 'data/result_fit.rda')
save(result_dat, file = 'data/results_df.rda')
