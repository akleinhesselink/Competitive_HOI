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

load('data/rs.rda')
load('data/rs_results.rda')
load('data/parms.rda')
load('data/eqs.rda')

nspp <- length(parms$r)

seeds <- do.call(rbind, rs_results)
rs <- do.call(rbind, rs)

lambda <- apply( seeds[ apply(rs, 1, sum) == 1 , ], 2, max)

y <- seeds/rs                 # calculate fecundity in all rs 
data <- data.frame(rs, y)
names(data) <- c(paste0('N', 1:nspp), paste0('Y', 1:nspp))

form1 <- paste('~ -1 + N1 + N2 + N3')
form2 <- paste(form1, '+ I(N1*N2) + I(N1*N3) + I(N2*N3)')

all_forms <- lapply( ls() [ ls() %>% str_detect('form') ] , function(x) eval(parse( text = x) ) )

# loop through species 
# prepare data for each species
# loop through models 
get_data <- function(x, spp){ 
  out <- x %>% select(starts_with('Y'))
  y <- out[, spp]
  dat <- as.matrix(x %>% select(starts_with('N')))
  dat[ , spp ] <- dat[ , spp] - 1
  out <- data.frame(y = y, dat)
  return(out[ complete.cases(out), ])
}

dat <- lapply( 1:nspp, get_data, x = data )
names( dat ) <- 1:nspp

dat <- do.call(rbind, dat)

dat <- dat %>% 
  tibble::rownames_to_column('id') %>% 
  separate( id, c('species', 'row'), sep = '\\.')

dat %>% filter(N3 == 0, N2 == 0 ) %>% ggplot( aes( x = N1, y = y)) + geom_point() + facet_wrap( ~ species)
dat %>% filter(N3 == 0, N1 == 0 ) %>% ggplot(aes(x = N2, y = y)) + geom_point() + facet_wrap( ~ species)
dat %>% filter(N2 == 0, N1 == 0 ) %>% ggplot( aes( x = N3, y = y)) + geom_point() + facet_wrap( ~ species)

fit_model <- function( dat, inits, model, frm, ...) { 
  
  fit <- optim( par = inits, fn = model, data = dat, form = frm, ...)
  pred <- model(fit$par, data = dat, form = frm, predict = T)
  return( list (fit = fit, data = data.frame( dat, pred = pred) ))
}

forms <- lapply( all_forms, as.formula)

dat <- 
  dat %>% mutate( num_sp = (N1 != 0) + (N2 != 0) + (N3 != 0) ) %>% 
  mutate( one_species = num_sp %in% c(0, 1)) %>%
  mutate( two_species = num_sp %in% c(0, 1, 2)) %>%
  mutate( three_species = num_sp %in% c(0:3)) 

n  <- 3
test <- dat %>% 
  filter( one_species) %>% 
  filter_( paste0( 'species == ', n ))

out <- fit_model(test , c(lambda[3], 0,0,0, -1), model = mod_bh, frm = forms[[1]], method = 'BFGS')
out$fit

out <- fit_model(test , c(lambda[3], 1,1,1, 1,1,1), model = mod_bh2, frm = forms[[1]], method = 'BFGS')
out$fit

out$data %>% 
  gather( comp, num , N1:N3 ) %>% 
  filter( num_sp == 0 | num > 0) %>% 
  ggplot( aes( x = num, y = y, color = comp)) + 
  geom_point() + 
  geom_line(aes( y = pred), linetype = 2)

# todo make this a function! 

sp <- 1
test <- dat %>% 
  filter( type ) %>% 
  filter_( paste0( 'species == ', sp ))

init <- seq(0, 1, by = 0.5)
  
inits <- expand.grid( c(1,1,1, rep( list(init), 3) )) 
inits <- inits %>% unique()
  
out1 <- fit_model(test , c(lambda[n], 1,1,1, -1), model = mod_bh, frm = forms[[1]], method = 'BFGS')
out1$fit

init <- seq(0, 1, by = 0.5)
inits <- expand.grid( c(1,1,1, rep( list(init), 3) )) 
inits <- inits %>% unique()

out2 <- fit_model(test , c(lambda[n], 1,1,1, 1,1,1), model = mod_bh2, frm = forms[[1]], method = 'BFGS')
out2$fit


init <- seq(0, 1, by = 0.5)
inits <- expand.grid( c(1,1,1,1,1,1,1,1,1, rep( list(init), 3) )) 
inits <- inits %>% unique()

tryfits <- list(NA)

for( i in 1:nrow(inits)) { 
  print(paste0('working on #', i))
  tryfits[[i]] <- try( fit_model(test , c(lambda[n], inits[i, ]), model = mod_bh2, frm = forms[[2]], method = 'BFGS', control = list( reltol = 1e-10)), silent = T)
  
  if(class(tryfits[[i]]) == 'try-error'){ 
    tryfits[[i]] <- NULL
  }
}

tryfits[[8]]$data

which.min( unlist( lapply ( tryfits, function(x) x$fit$value) ))

tryfits[[1]]$fit

out1$data %>% ggplot( aes( x = y, y = pred)) + geom_point()
out3$data %>% ggplot( aes( x = y, y = pred)) + geom_point()


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
