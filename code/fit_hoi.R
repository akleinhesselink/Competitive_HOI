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

load('data/experiments.rda')
load('data/out.rda')
load('data/parms.rda')
load('data/seedling_mass.rda')
load('data/conversion.rda')

fecundity <- phenology <- data.frame( matrix( NA, nrow = nrow(experiments), ncol = 3))

get_resource_use <- function(State, f, parms){ 
  State[3:5]*f(State[2], parms$r, parms$K)
}

use <- lapply( out, function(x) t(apply( x, 1, get_resource_use, f, parms)))
experiments <- split(experiments,1:nrow(experiments))
get_per_capita_use <- function(x, y) sweep( x, 2, as.numeric(y), '/')
per_capita_use <- mapply(use, experiments, FUN = get_per_capita_use, SIMPLIFY = F)
phenology <- lapply( out, function( x ) apply(x[, c(3:5)], 2, find_phenology))

##

experiments <- do.call(rbind, experiments)
max_biomass <- do.call( rbind, lapply( out, function( x ) apply( x[, c(3:5)], 2, max)))
fecundity <- (max_biomass*conversion)/seedling_mass
phenology <- do.call(rbind, phenology )

rm(out)

lambda <- diag( as.matrix( fecundity [ apply(experiments, 1, sum) == 1, ]  )) # calculate lambdas 
y <- fecundity/experiments                                                    # calculate fecundity in all experiments 
data <- data.frame(experiments, y, phenology/10)
names(data ) <- c('N1', 'N2', 'N3', 'Y1', 'Y2', 'Y3', paste0('PH', c(1:3)))

form1 <- paste('~ -1 + N1 + N2 + N3')
form2 <- paste(form1, '+ I(N1*N2) + I(N1*N3) + I(N2*N3)')

all_forms <- c(form1, form2)

gg_intra <- gg_inter <- gg_phenology <- list()
results2 <- results <-  list()
i = 1

# fit annual plant model parameters --------------------------------------- # 
for(i in 1:3){ 
  # loop through each species and fit annual plant model parameters ------------------------------------------- # 
  data1 <- data
  data1[, i] <- data1[, i] - 1            # remove focal from competive neighborhood 
  data1$y  <- data1[, i + 3]              # focal per capita seed production
  data1$data <- as.matrix(data1[,c(1:3)])
  data2 <- data1[is.finite(data1$y),]  
  
  data2 <- data2[ (data2[, i] == 0 & (data2[, c(1:3)[-i][1]] == 0 & data2[, c(1:3)[-i][2]] == 0 )) | (data2[,i] > 0 & sum(data2[,c(1:3)[-i]]) > 1), ]
  
  terms <- lapply( all_forms, function(x) length(str_split(x, pattern = '\\+')[[1]] ))
  forms <- lapply(all_forms, as.formula)
  
  pars <- lapply( terms, function(x, ...) {c( ..., rep(1, x - 1), -1 )} , Z = lambda[i] )  
  
  out <- mapply(par = pars, form = forms, FUN = optim, MoreArgs = list(data = data2, fn = mod_bh, method = 'BFGS', control = list( maxit = 1e9, reltol = 1e-15)), SIMPLIFY = F)
  results[[i]] <- out 
  
  pars <- lapply(out, function(x) x$par)
  predictions <- mapply(par = pars, form = forms, FUN = mod_bh, MoreArgs = list(data = data1, predict = T))
  predictions <- data.frame(predictions)  
  names(predictions) <- paste0(paste0('pY', i , '.model'), 1:length(forms))
  data <- cbind(data, predictions)
  
} 


data <- 
  data %>% 
  gather(prediction, value,  starts_with('p', ignore.case = F)) %>% 
  gather( species, observed, Y1:Y3) %>% 
  separate( prediction, c('pred_species', 'model'), sep = '\\.') %>% 
  filter(str_extract(pred_species, 'Y\\d+') == str_extract(species, 'Y\\d+')) %>% 
  select( - pred_species) %>% 
  spread( model, value ) %>% 
  filter( is.finite(observed)) %>% 
  mutate( diff1 = (model1 - observed)/observed, diff2 = (model2 - observed)/observed) %>% 
  gather( model, deviation, diff1, diff2)


ggplot( data %>% filter( N1 == 20, N2 > 15, N3 > 15), aes( x = N2, y = N3, fill = deviation )) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  facet_wrap(~ model)

ggplot( data %>% filter( N2 == 25, N1 > 15, N3 > 15), aes( x = N1, y = N3, fill = deviation )) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  facet_wrap(~ model)

ggplot( data %>% filter( N3 == 30, N1 > 15, N2 > 15), aes( x = N1, y = N2, fill = deviation )) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  facet_wrap(~ model)


basic_plot <- function(data, focal, hold_at = 4, c1range, c2range) { 
  N_focal <- paste0('N', focal)
  Y_focal <- paste0('Y', focal)
  comp <- paste0('N', c(1:3)[-focal])
  
  data <- 
    data %>% 
    filter( species == Y_focal) %>% 
    filter_(.dots = paste0(N_focal, ' == ', hold_at)) %>% 
    filter_(.dots = paste0(comp[1], ' %in% ', c1range )) %>% 
    filter_(.dots = paste0(comp[2], ' %in% ', c2range )) %>% 
    select( N1:N3, species, observed, model1, model2) %>% 
    gather(type, predicted, model1:model2)
  
  data[, comp[2]] <- as.factor(data[, comp[2]])
  
  ymx <- max( c(data$observed, data$predicted))
  
  data %>% 
    ggplot(aes_string( x = comp[1], y = 'observed', color = comp[2])) +
    geom_point() +
    geom_line( aes(y = predicted, linetype = type)) + 
    facet_wrap( ~ type ) 
}

data %>% basic_plot(focal = 1, hold_at = 20, c1range = '25:35', c2range = '30:40' )
data %>% basic_plot(focal = 2, hold_at = 30, c1range = '15:25', c2range = '30:40' )
data %>% basic_plot(focal = 3, hold_at = 35, c1range = '15:25', c2range = '25:35' )

# compare error of fits 
fits <- data.frame( do.call( rbind, ( lapply( results, function(x) unlist(lapply(x, function(x) x$value))))))
fits$species <- paste('species', 1:3)
fits <- fits %>% gather( model, 'MSE', X1:X2, -species)
fits$model_label <- factor(fits$model, labels = c('m1', 'm2'))
fits <- fits %>% group_by( species ) %>% mutate( MSE_rel = MSE/max(MSE))

ggplot( fits %>% filter( model == 'X2' ), aes( x = species, y = MSE_rel)) + 
  geom_bar(stat = 'identity')

fit_plot <- ggplot( fits, aes( x = model_label, y = MSE) ) + 
  geom_bar(stat= 'identity') + 
  facet_grid(species ~. , scales = 'free') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5 ), axis.title.x = element_blank()) 

fit_plot

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
