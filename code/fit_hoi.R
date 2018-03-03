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


fill_in_inits <- function(x) {  
  if( length(x) > 1) { 
    seq(x[1], x[2], x[3]) 
  }else if( length(x) == 1) { 
    x
  }
}


fit_inits <- function(inits, model, dat, frm, lower_b, ... ){ 
  
  inits <- expand.grid( inits ) 
  j <- nrow(inits)
  fits <- list(NA)
  
  inits <- split( inits, 1:nrow(inits))
  
  fits <- mclapply( inits,  function(x){ 
    out <- try( fit_model(dat , x, model = model, frm = frm, method = 'L-BFGS-B', lower=lower_b), silent = T)
    if(class(out) == 'try-error'){ 
      out <- NULL
    }
    return(out)
  }, ... )
  
  fits <- fits[ !unlist(lapply( fits, is.null) ) ]  
  
  fits <- unlist( fits[ which.min( lapply( fits, function(x) x$fit$value) ) ], recursive = F)
  names(fits) <- c('fit', 'data')
  
  return( fits )
}

model_fits <- list(NA)
n = 1
for( n in 1:nspp) { 
  
  temp_dat <- dat %>% 
    filter(one_species) %>% 
    filter_( paste0( 'species == ', n ))

  inits <- list(lambda[n], 1, 1, 1, c(-1.1,-0.1,0.5))
  lower_b <- c(1, rep(0,length(inits) -2 ), -3)
  res1 <- fit_inits(lapply( inits, fill_in_inits), mod_bh, temp_dat, forms[[1]], lower_b, mc.cores = 4)

  inits <- list(lambda[n], c(0.01, 1.01, 0.5), c(0.01, 1.01, 0.5), c(0.01, 1.01, 0.5))
  lower_b <- c(1, rep( 0, length(inits) - 1))
  res2 <- fit_inits(lapply( inits, fill_in_inits), mod_bh2, temp_dat, forms[[1]], lower_b, mc.cores = 4)
  
  inits <- list(lambda[n], 1, 1, 1, 1, 1, 1, c(-1.1, -0.1, 0.5))
  lower_b <- c(1, rep(0,length(inits) -2 ), -3)
  res3 <- fit_inits(lapply( inits, fill_in_inits), mod_bh, temp_dat, forms[[2]], lower_b, mc.cores = 4)

  inits <- list(lambda[n], 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  lower_b <- c(1, rep(0,length(inits) - 1))
  res4 <- fit_inits(lapply( inits, fill_in_inits), mod_bh2, temp_dat, forms[[2]], lower_b, mc.cores = 4)

  model_fits[[n]] <- list(m1 = res1, m1HOI = res2, m2 = res3, m2HOI = res4)

}


names(model_fits) <- paste0('N', c(1:nspp))


save(model_fits, file = 'data/model_fits.rda')


res1$data %>% head()


