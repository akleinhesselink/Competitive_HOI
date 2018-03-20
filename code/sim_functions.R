library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)
library(gridExtra)
library(scales)

# annual plant model formula 

form1 <- as.formula('~ -1 + N1 + N2 + N3')
formHOI <- as.formula('~ -1 + N1 + N2 + N3 + I(N1*N2) + I(N1*N3) + I(N2*N3)')


f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, r, K, q, m) { B*(q*f(R, r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K, p, epsilon) { p[u] - epsilon*R - sum(B*f(R,r, K)) } # resource (water)


TO_fun <- function(r, TO_pars) {
  # trade-off between uptake at low resource availability and uptake during high resource availability
  with(TO_pars, { 
    alpha + (beta - alpha)*((r - gamma)/(delta - gamma))^cc
  })
}

Rstar <- function(r, K, m, q) { m*K/(q*r-m) }      # resource required for growth to balance loss
find_phenology <- function(B) { min(which(B == 0)) } # find the flowering time

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, r, K, q, m)
    dR <- dRdu(u, B, R, r, K, p, epsilon)
    return( list( c(dR, dB))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K/(q*r-m) < 0.000001) # logical vector of species to terminate
    State[2:length(State)][ terminate ] <- 0
    return(State)
  })
}

run_experiment <- function(seedlings, parms) { 
  seedlings <- as.numeric(seedlings)
  State <- c(parms$soil_m, seedlings*parms$seedling_mass)
  out <- ode(y=State, times = seq( 1, parms$times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root)
  max_biomass <- apply( out[, 3:5 ], 2,  max )
  
  ( max_biomass*parms$conversion )/parms$seedling_mass
}

run_multi_gen <- function(seedlings, t, parms, tol){ 
  
  out <- list(NA)
  population <- data.frame(N1 = rep(NA, t), N2 = rep(NA, t), N3 = rep(NA, t))
  population[1, ] <- seedlings
  
  i <- 1
  pop_diff <- rep(1, 3)
  while( i < t & any(pop_diff > tol)){
    State <- as.numeric(c(parms$soil_m, population[i, ]*parms$seedling_mass))
    out[[i]] <- ode(y=State, times = seq(1, parms$times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
    max_biomass  <- apply( out[[i]][, c(3:5)], 2, max)
    population[i+1,] <- (max_biomass*parms$conversion)/parms$seedling_mass
    pop_diff <- abs( population[i+1,] - population[i, ])
    i <- i + 1 
  }
  
  return( population) 
}

plot_timeseries <- function(x, parms, ... ){ 
  par(mfrow = c(2,1))
  nspp <- ncol(x) - 2
  with(parms, {
    flowering_times <- x[apply( x[, c(3:(2+nspp)), drop = F ], 2, find_phenology), 1]
    plot( x[,1], x[,2], type = 'l', ylim = c(0,max(x[,2])), ylab = 'Soil moisture', xlab = 'day')
    abline( h = Rstar(r=r, K = K, m = m, q = q), lty = 2, ...)
    matplot(x = x[,1], x[, c(3:(2+nspp))], type = 'l', xlab = 'day', ylab = 'biomass', lty = 1, ...)
  })
  legend(5, -0.1, legend = paste('species', 1:nspp), col = my_colors[1:nspp], cex = 1, xpd = T, lty = 1, bg = NA, box.col = NA, yjust = 0)
}

plot_transpiration <- function(parms, my_colors ){
  # plot transpiration rate curve 
  par(mfrow =c(1,1))  
  nspp <- length(parms$r)
  with(parms, {
    plot(R/500, f(R, r=r[1], K = K[1]), type = 'l', ylab = 'Transpiration rate (mm water per g per day)', xlab = 'Soil moisture content (%)')
    lapply(1:nspp, function(x,...){ points(R/500, f(R, r=r[x+1], K = K[x+1]), type ='l', col = my_colors[x+1])})
    legend( 'bottom', legend = paste('species', 1:nspp), col = my_colors[1:nspp], cex = 1, xpd = T,  lty = 1, bg = NA, box.col = NA)
  })
}

plot_growth_rate <- function(parms, my_colors){
  # plot relative growth rate curves 
  par(mfrow = c(1,1))
  nspp <- length(parms$r)
  with(parms, {
    plot(R, q*f(R, r=r[1], K = K[1]) - m , type = 'l', ylab = 'Growth rate (g per g per day)', xlab = 'Soil moisture (mm^3 of water per cm^3)')
    lapply(1:nspp, function(x,...){ points(R, q*f(R, r=r[x+1], K = K[x+1]) - m, type = 'l', col = my_colors[x+1])})
    abline(h = 0, lty = 2)
    legend( 'bottom', legend = paste('species', 1:nspp), col = my_colors[1:nspp], cex = 1, xpd = T,  lty = 1, bg = NA, box.col = NA)
  })
}

plot_Rstar <- function(parms, my_colors){ 
  par(mfrow = c(1,1))
  nspp <- length(parms$r)
  with(parms, {
    plot(Rstar(r,K,m,q),q*f(1e8, r, K) - m,  type = 'l', xlab = 'R*', ylab = 'Max growth rate')
    lapply(1:nspp, function(x, ... ){ points(Rstar(r[x],K[x],m,q),q*f(1e8, r[x], K[x]) - m ,  col = my_colors[x], pch = 20, cex = 2) }) 
    legend( 'topleft', legend = paste('species', 1:nspp), col = my_colors[1:nspp], cex = 1.2, xpd = T,  pch = 20, bg = NA, box.col = NA)
  })
}


seed_production <- function(x, init_size, R_state, parms, func, TO_fun, TO_pars, ...){ 
  x <- matrix(x, 1, length(x))
  State <- c(R_state, x*init_size)
  parms$K <- TO_fun(parms$r, TO_pars)
  out <- ode(y=State, func = func, parms = parms, ...)
  peak_biomass <- apply( out[, c(3:ncol(out)), drop = F], 2, max)
  seeds_out <- peak_biomass*conversion
  return(seeds_out)
}

fitness <- function(x, fun, ...){ 
  log(fun(x, ...)/x) 
}

find_N_hat <- function(my_range = c(1e-7, 100), fun, ... ){ 
  obj_fun <- function(x, ...){ abs(fitness(x, fun, ... )) }
  optimize(obj_fun, interval = my_range, ... )
}

find_N_hat2 <- function(my_range = c(1e-7, 100), fun, ... ){ 
  uniroot(fitness, interval = my_range, fun = fun, ... )
}

ann_plant_mod <- function(x, form, pars) { 
  pars <- unlist(pars)
  mm <- model.matrix(as.formula( form ), x )
  nt <- ncol(mm)
  nsp <- ncol( x %>% select(starts_with('N')))
  
  lambda <- head( pars, nsp)   
  tau <- tail(pars, nsp)
  alphas <- pars[ -which(pars == lambda | pars == tau) ]

  alphas <- matrix( alphas, nsp, nt)

  y <- lambda*((1 + rowSums( sweep(alphas, 2, mm, '*') ))^(tau))    
  
  return(y)
}

fit_ann_plant <- function(data, model, form, focal = 1, my_inits = NULL, ... ){ 
  
  temp <- 
    data %>% 
    filter_( .dots = paste0( 'focal == "F', focal, '"')) 
  
  mm <- model.matrix(form, temp)
  
  npar <- ncol(mm)

  if( is.null(my_inits) ){ 
    par <- c( max(temp$fecundity), rep(1, npar), -1)
  }else{ 
    par <- my_inits 
  }

  optim(par = par, fn = model, y = temp$fecundity, mm = mm, ... )
}

predict_fit <- function( dat, model, pars, form, foc = 1 ){ 
  mod_name <- deparse(substitute(model))
  form_name <- deparse(substitute(form))
  dat <- 
    dat %>% 
    filter( focal == paste0('F', foc) )

  dat$y <- NA
  mm <- model.matrix( form, dat)
  
  dat$pred <- as.numeric(model(pars = pars, dat$y, mm, predict = T))

  dat <- dat %>% 
    ungroup() %>% 
    select(id, focal, pred) %>% 
    filter( !is.na(pred)) %>% 
    distinct() %>% 
    arrange( as.numeric(id)) %>% 
    mutate( focal_predicted = paste0('pred.', mod_name, '.', focal)) %>% 
    spread( focal_predicted, pred )
  
  dat$form <- form_name
  
  dat
  
}


make_monoculture <- function(experiments) { 
  
  if( !identical( names(experiments), c('N1', 'N2', 'N3', 'id')) ) { 
    stop( 'incorrect format for experiments') 
  }
  
  experiments <- 
    experiments %>%   
    arrange(as.numeric(id)) %>%
    filter( (N1 == 0 & N2 == 0) | (N3 == 0 & N2 == 0 ) | (N1 == 0 & N3 == 0 ) )
  
  return(experiments)  
}

make_biculture <- function(experiments) { 
  
  if( !identical( names(experiments), c('N1', 'N2', 'N3', 'id')) ) { 
    stop( 'incorrect format for experiments') 
  }
  
  experiments <- 
    experiments %>%   
    arrange(as.numeric(id)) %>%
    filter( (N1 == 0) | (N2 == 0 ) | (N3 == 0 ) )
  
  return(experiments)  
}

make_experiments <- function(maxdens = 10, base = 2, nspp) { 
  # -------- simulate gradient  -------------------------- # 
  experiments <- expand.grid(N1 = c(0, c(base^c(0:maxdens))), N2 = c(0, base^c(0:maxdens)), N3 = c(0, base^c(0:maxdens)))
  experiments$id <- row.names(experiments)
  return(experiments)
}

add_focal <- function(X, focal, focal_lab) { 
  out <- X
  out[, grep('N', names(out)) ] <- data.frame( t(apply( out[, grep('N', names(out)) ], 1, function(x) x  + focal)))
  out$focal <- focal_lab 
  out
}

basic_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

plot_two_sp <- function( data, focal = 'F1', C1 = 'N1', C2 = 'N2', C3 = 'N3', C2_dens = c(0, 16, 64, 256, 1024) ) { 
  
  data <- data[ data$focal == focal, ] 
  
  data$C1 <- as.numeric(unlist(data[, C1]))
  data$C2 <- as.numeric(unlist(data[, C2]))
  data$C3 <- as.numeric(unlist(data[, C3]))
  
  data %>%
    filter( comp_n == 0 | ( C3 == 0 )) %>%
    select(-C3) %>%
    filter( C2 %in% C2_dens ) %>%  
    mutate( C2 = as.factor(C2)) %>%
    ggplot( aes( y = fecundity, x = C1, color = C2)) +
    geom_point() +
    scale_y_continuous(name = paste0('N', str_extract(focal, '\\d'), ' fecundity'), trans = 'log', breaks = basic_breaks()) +
    scale_x_continuous(name = paste(C1, 'density')) + 
    scale_color_discrete(name = paste(C2, 'density'))
}

mod_bh <- function(pars, y, mm, predict = FALSE){ 
  
  if( !length(pars) == ncol(mm) + 2 ){ stop('wrong number of parameters supplied!')}
  
  lambda <- pars[1]
  tau <- pars[length(pars)]
  alphas <- pars[2:(length(pars) - 1)]
  
  mu <- lambda*(1 + mm%*%alphas)^tau
  
  Error <- sum( (mu - y )^2 )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(Error) 
  }
}


mod_bh2 <- function( pars, y, mm, predict = F){ 
  
  if(!length(pars) == ncol(mm) + 1 ){ stop('wrong number of parameters supplied!')}
  
  alphas <- pars[ 2:length(pars) ]
  lambda <- pars[1]
  
  mu <- NA
  for( i in 1:nrow(mm)){ 
    mu[i] <- lambda/( 1 + sum(  mm[i, ]^alphas  ) )
  }
  
  Error <- sum( (mu - y)^2 ) 
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(Error) 
  }
  
}

mod_bh_ll <- function(pars, y, mm, sd = 1, predict = FALSE){ 
  
  if( !length(pars) == ncol(mm) + 2 ){ stop('wrong number of parameters supplied!')}
  
  lambda <- pars[1]
  tau <- pars[length(pars)]
  alphas <- pars[2:(length(pars) - 1)]
  
  mu <- lambda*(1 + mm%*%alphas)^tau
  
  neg_ll <- sum( - dnorm( mu, y, sd, log = T) )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(neg_ll) 
  }
}

mod_bh2_ll <- function(pars, y, mm, sd = 1, predict = FALSE){ 
  
  if( !length(pars) == ncol(mm) + 1 ){ stop('wrong number of parameters supplied!')}
  
  lambda <- pars[1]
  alphas <- pars[2:(length(pars))]
  
  mu <- NA
  for( i in 1:nrow(mm)){ 
    mu[i] <- lambda/( 1 + sum(  mm[i, ]^alphas  ) )
  }
  
  neg_ll <- sum( - dnorm( mu, y, sd, log = T) )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(neg_ll) 
  }
}


fit_2_converge <- function(n_seq, start_sd, min_sd,  my_inits, model, ...){ 
  res <- list()
  inits <- list()
  inits[[1]] <- my_inits
  
  sd_grad <- rev( seq(min_sd, start_sd, length.out = n_seq))
  sd <- sd_grad[1]
  res[[1]]  <- fit_ann_plant(sd = sd, my_inits = my_inits, model = model, ...)
  
  for( i in 2:n_seq) { 
    sd <- sd_grad[i]
    inits[[i]] <- res[[i-1]]$par
    res[[i]] <- fit_ann_plant( my_inits = inits[[i]], sd = sd, model = model, ... )
  }
  converged <- unlist( lapply( res, function(x) x$convergence == 0 ))
  sd_grad <- sd_grad[which(converged)]
  res <- res[which(converged)]
  
  return( list( sd_grad = sd_grad, res = res) )
}  

fit_both_mods <- function( focal = 1, form1, inits1, lower1, model, ... ){  
  mod_name <- deparse(substitute(model))
  
  fits <- fit_2_converge(focal = focal, form = form1, my_inits = inits1, lower = lower1, model, ... )
  
  best_fit <- which.min( unlist( lapply( fits$res, function(x)  x$value ) ) )
  fit1 <- fits$res[[best_fit]]
  lambda <- fit1$par[1]
  
  if(str_detect(mod_name, '2')) { 
    alpha  <- 2:length(inits1)
    alpha  <- fit1$par[ 2:length(fit1$par) ]
    
    HOI_inits <- c(lambda, alpha, inits1[-1])
    lower  <- c(lower1[1], lower1[-1])
    
  }else{ 
    alpha  <- fit1$par[-c(1, length(fit1$par)) ]
    tau    <- fit1$par[length(fit1$par)]  
    HOI_inits <- c(lambda, rep(0, 2*length(alpha)), tau)
    lower  <- rep( 0, length(HOI_inits))
    lower[length(lower)] <- -1 
  }
  
  lower[1] <- lower1[1]
  
  fits <- fit_2_converge(focal = focal, form = formHOI, my_inits = HOI_inits, lower = lower, model, ... )
  
  best_fit <- which.min( unlist( lapply( fits$res, function(x)  x$value ) ) )
  fitHOI <- fits$res[[best_fit]]
  
  return( list(fit1 = fit1, fitHOI = fitHOI))
  
}

