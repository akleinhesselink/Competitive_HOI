
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

mod_bh <- function(pars, data, form, predict = FALSE){ 
  
  if(length(pars) < 2){ stop('not enough parameters supplied!')}
  
  mm <- model.matrix(form, data)
  
  mu <- pars[1]*(1 + mm%*%pars[c(2:(length(pars) - 1))])^pars[length(pars)]
  
  Error <- sum( (log(mu) - log(data$y))^2 )

  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(Error) 
  }
}

mod_bh3 <- function(pars, data, form, predict = FALSE){ 
  
  if(length(pars) < 2){ stop('not enough parameters supplied!')}
  
  mm <- model.matrix(form, data)
  
  mu <- pars[1]*(pars[2] + mm%*%pars[c(3:(length(pars) - 1))])^pars[length(pars)]
  
  Error <- sum( (log(mu) - log(data$y))^2 )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(Error) 
  }
}





mod_bh2 <- function( pars, data, form, predict = F, lg = F){ 
  
  if(length(pars) < 2){ stop('not enough parameters supplied!')}
  
  mm <- model.matrix(form, data)
  lambda <- pars[1]
  pars <- pars[-1]
  
  taus <- pars[ 1:ncol(mm) ]

  mu <- NA
  for( i in 1:nrow(mm)){ 
    mu[i] <- lambda/( 1 + sum(  mm[i, ]^taus  ) )
  }
  
  Error <- sum( (log(mu) - log(data$y))^2 ) 

  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(Error) 
  }
  
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
  alphas <- pars[ rev(rev(seq_along(pars))[-c(1:nsp)])[-c(1:nsp)] ]
  alphas <- matrix( alphas, nsp, nt, byrow = T)
  
  y <- lambda*((1 + mm%*%alphas)^(tau))    
  
  return(y)
}

fit_ann_plant <- function(data,  model, focal = 1, my_inits = NULL, ... ){ 
  
  temp <- 
    data %>% 
    filter_( .dots = paste0( 'focal == "F', focal, '"')) %>% 
    spread( competitor, density, fill = 0)
  
  nspp <- ncol(temp %>% select(starts_with('N')))
  
  temp$y <- temp$fecundity
  if( is.null(my_inits)){ 
    par <- c(max(temp$fecundity), rep(1, nspp), -1)
  }else{ 
    par <- my_inits 
  }
  optim(par = par, fn = model, data = temp, ... )
}

predict_fit <- function( dat, model, pars, form, foc = 1 ){ 
  mod_name <- deparse(substitute(model))
  
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
    mutate( focal_predicted = paste0('pred.', mod_name, '.', focal)) %>% 
    spread( focal_predicted, pred )
}


make_monoculture <- function(results) { 
  
  if( !identical( names(results), c('id', 'N1', 'N2', 'N3', 'F1', 'F2', 'F3')) ) { 
    stop( 'incorrect format for results') 
  }
  
  results <- 
    results %>%   
    arrange(as.numeric(id)) %>%
    filter( (N1 == 0 & N2 == 0) | (N3 == 0 & N2 == 0 ) | (N1 == 0 & N3 == 0 ) ) %>% 
    mutate( lambda =  ifelse(N1 == 0 & N2 == 0 & N3 == 0 , T, F)) %>% 
    gather( competitor, density, N1:N3) %>% 
    filter( lambda | density > 0 ) %>% 
    gather( focal, fecundity, F1:F3)  
  
  results$competitor_label <- paste0( 'competitor\n', results$competitor) 
  results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))
  
  return(results)  
}

fit_2_converge <- function(results, model, ... ){ 
  fits <- lapply(1:nspp, function(x, ...) fit_ann_plant(focal = x, ... ), data = results, model = model, ... )
  converged <- lapply( fits, function(x) x$convergence) == 0
  
  while ( ! all(converged) ) { 
    fitpars <- lapply( fits, function(x) x$par )
    fits[!converged] <- mapply(x = which(!converged), y = fitpars[!converged], FUN = function(x, y, ... ) fit_ann_plant(focal = x, my_inits = y, model = model, ... ), MoreArgs = list( data = results, ...), SIMPLIFY = F)
    converged <- lapply( fits, function(x) x$convergence) == 0
  }
  
  return(fits)
}

form1 <- as.formula('~ -1 + N1 + N2 + N3')
formHOI <- as.formula('~ -1 + N1 + N2 + N3 + I(N1*N2) + I(N1*N3) + I(N2*N3)')


make_experiments <- function(maxdens = 10, base = 2, nspp) { 
  # -------- simulate gradient  -------------------------- # 
  experiments <- expand.grid(N1 = c(0, c(base^c(0:maxdens))), N2 = c(0, base^c(0:maxdens)), N3 = c(0, base^c(0:maxdens)))
  experiments$id <- row.names(experiments)
  return(experiments)
}


