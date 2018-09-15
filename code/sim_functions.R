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

plot_timeseries <- function(out, sp_labs = c('Resource', '1','2'), mytheme, fname = 'figures/example_timeseries.png'){ 
  library(gridExtra)
  
  temp <- 
    data.frame(out) %>% 
    gather( var, val, starts_with('X')) %>% 
    mutate( species = var ) %>% 
    mutate( species = factor( species, labels = sp_labs )) %>% 
    filter( time < 150 ) 
    
  resource_plot <- 
    ggplot( temp %>% filter( species == 'Resource'), aes( x = time, y = val )) + 
    geom_line() + 
    ylab( 'Resource') + 
    mytheme + 
    theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.text = element_blank()) 
  
  biomass_plot  <- ggplot( temp %>% filter( species != 'Resource'), aes( x = time, y = val, color = species)) + 
    geom_line() + 
    ylab( 'Biomass') + 
    xlab( 'Day of Year') + 
    scale_color_manual(values = my_colors, guide = F) + 
    mytheme +
    theme(legend.position = c(0.85, 0.5), 
          legend.background = element_rect(colour = 1, size = 0.2), 
          axis.text = element_blank()) 
    

  return(list(resource_plot, biomass_plot))
}

plot_transpiration <- function(parms, my_colors ){
  # plot transpiration rate curve 
  par(mfrow =c(1,1))  
  nspp <- length(parms$r)
  with(parms, {
    plot(R/500, f(R, r=r[1], K = K[1]), type = 'l', ylab = 'Resource uptake rate (units of R per g per day)', xlab = 'Soil resource content (%)')
    lapply(1:nspp, function(x,...){ points(R/500, f(R, r=r[x+1], K = K[x+1]), type ='l', col = my_colors[x+1])})
    legend( 'bottom', legend = paste('species', 1:nspp), col = my_colors[1:nspp], cex = 1, xpd = T,  lty = 1, bg = NA, box.col = NA)
  })
}

plot_resource_uptake <- function(parms, R = 0:500, spec_labs = c('1','2')){ 
  
  curves <- data.frame(R = R,  mapply(x = as.list(parms$r), y = as.list(parms$K), FUN = function(x, y) { f(R = R, x, y) }) )
  
  curves$X1[curves$R < Rstar(parms$r[1], parms$K[1], parms$m, parms$q) ] <- 0
  curves$X2[curves$R < Rstar(parms$r[2], parms$K[2], parms$m, parms$q) ] <- 0
  curves$X3[curves$R < Rstar(parms$r[3], parms$K[3], parms$m, parms$q) ] <- 0
  
  curves <- 
    curves %>% 
    gather( species, uptake, starts_with('X')) %>% 
    filter( uptake > 0 ) %>% 
    mutate( species = factor(species, labels = spec_labs))
  

  curves %>% 
    ggplot(aes( x = R, y = uptake, color = species )) + 
    geom_line() +
    my_theme +
    ylab( 'Resource uptake rate \n per g per day') + 
    xlab( 'Resource') + 
    scale_color_manual(values = my_colors, 'Species') + 
    theme(legend.position = c(0.75, 0.25), legend.background = element_rect(color = 1, size = 0.25))
  
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
    axisTicks(range(x, na.rm = TRUE), log = F, n = n)
  }
}


plot_two_sp <- function( data, focal = 'F1', C1 = 'N1', C2 = 'N2', C3 = 'N3', 
                         C2_dens = NULL ) { 
  
  data <- data[ data$focal == focal, ] 
  
  data$C1 <- as.numeric(unlist(data[, C1]))
  data$C2 <- as.numeric(unlist(data[, C2]))
  data$C3 <- as.numeric(unlist(data[, C3]))
  
  if(is.null(C2_dens)){
    C2_factor <- as.numeric( factor(data$C2)  )
    C2_dens <- floor( seq(min(C2_factor), max(C2_factor), length.out = 4))
  } 
  
  data %>%
    filter( comp_n == 0 | ( C3 == 0 )) %>%
    select(-C3) %>%
    mutate( C2 = as.factor(C2)) %>%
    filter( as.numeric(C2) %in% C2_dens ) %>%  
    ggplot( aes( y = fecundity, x = C1, color = C2)) +
    geom_point() +
    scale_y_continuous(name = paste0('N', str_extract(focal, '\\d'), ' fecundity'), breaks = basic_breaks()) +
    scale_x_continuous(name = paste(C1, 'density')) + 
    scale_color_discrete(name = paste(C2, 'density'))
}

mod_bh_ll <- function(pars, y, mm, sd = 1, predict = FALSE){ 
  
  if( !length(pars) == ncol(mm) + 2 ){ stop('wrong number of parameters supplied!')}
  
  lambda <- pars[1]
  tau <- pars[2]
  alphas <- pars[3:(length(pars))]
  
  mu <- lambda*(1 + mm%*%alphas)^tau
  
  neg_ll <- sum( - dnorm( log(mu), log(y), sd, log = T) )
  
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
  
  neg_ll <- sum( - dnorm( log(mu), log(y), sd, log = T) )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(neg_ll) 
  }
}

mod_bh3_ll <- function(pars, y, mm, sd = 1, predict = FALSE){ 
  
  if( !length(pars) == ncol(mm) + 1 ){ stop('wrong number of parameters supplied!')}
  
  lambda <- pars[1]
  alphas <- pars[2:(length(pars))]
  
  mu <- NA
  for( i in 1:nrow(mm)){ 
    mu[i] <- lambda*exp( - sum(mm[i, ]*alphas) )
  }
  
  neg_ll <- sum( - dnorm( log(mu), log(y), sd, log = T) )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(neg_ll) 
  }
}

mod_bh4_ll <- function(pars, y, mm, sd = 1, predict = FALSE){ 
  
  if( !length(pars) == ncol(mm)*2 + 1 ){ stop('wrong number of parameters supplied!')}
  
  lambda <- pars[1]
  tau <- pars[2:(ncol(mm) + 1)]
  alphas <- tail(pars, ncol(mm))
  
  mu <- NA
  for( i in 1:nrow(mm)){ 
    mu[i] <- lambda/( 1 + alphas%*%c(mm[i,]^tau))
  }
  
  neg_ll <- sum( - dnorm( log(mu), log(y), sd, log = T) )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(neg_ll) 
  }
}


compare_parameters <- function(original_file, fitted){ 
  
  original <- readRDS(original_file)

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
  
}


plot_all_fits <- function(all_fits){
  
  comp_spp <- names(all_fits) [ grep('^N\\d+', names(all_fits)) ] 
  focal_spp <- unique(all_fits$focal) 
  nspp <- length(focal_spp)
  
  fit_plot <- list()
  
  for( focal in 1:nspp){ 
    
    c1 <- sort(c(comp_spp[focal], comp_spp[-focal])[-nspp])
    c2 <- comp_spp[-focal]
    
    combos <- 
      expand.grid( focal = focal_spp[focal], 
                   c1 = c1, 
                   c2 = c2) %>% 
      filter( as.character(c1) != as.character(c2)) %>% 
      arrange(c1, c2)
    
    for(i in 1:nrow(combos)){ 
      combos$c3[i] <- comp_spp[ ! comp_spp %in% unlist((combos %>% select( c1, c2))[i, ] ) ]
    }
    
    combos <- as.matrix(combos)
    
    p <- list()
    for( i in 1:nrow(combos)){ 
      p[[i]] <- plot_two_sp(all_fits, 
                            focal = combos[i, 1], 
                            C1 = combos[i, 2], 
                            C2 = combos[i, 3], 
                            C3 = combos[i, 4]) +
        geom_line(aes( y = pred_fecundity, linetype = form)) + 
        scale_linetype_manual(values = c(2,3), guide = F) + 
        theme(legend.position = c(0.85, 0.7))
      
      if(i == nrow(combos)){ 
        p[[i]] <- plot_two_sp(all_fits, 
                    focal = combos[i, 1], 
                    C1 = combos[i, 2], 
                    C2 = combos[i, 3], 
                    C3 = combos[i, 4]) +
          geom_line(aes( y = pred_fecundity, linetype = form)) + 
          scale_linetype_manual(values = c(2,3)) + 
          theme(legend.position = c(0.85, 0.7))
      }
      
      
    }
    
    fit_plot[[focal]] <- do.call( function(...) grid.arrange (..., ncol = nspp), p )
    
  } 
  
  fit_plot
}


fit_model <- function(dat, form = form1, mod_name = "mod_bh_ll", start_sd = 3, min_sd = 1, max_refit = 5){
  
  model <- eval(parse(text = mod_name))
  
  sd_grad <- rev(seq(min_sd, start_sd, length.out = max_refit))
  form_name <- deparse(substitute(form))

  if( str_detect(form_name, 'HOI') ){
    if( str_detect(mod_name, '4')){ 
      # HOI TYPE 4 
      init <- c(20, rep(1,6), rep(0, 6))
      lower <- c(1, rep(0,6), rep(0, 6))
      upper <- c(1e2, rep(2, 12))
    }else if( str_detect(mod_name, '3')){ 
      # HOI TYPE 3 
      init <- c(20, 1, 1, 1, 0, 0, 0)
      lower <- c(1, 0, 0, 0, -0.01, -0.01, -0.01)
      upper <- c(1e2, 2, 2, 2, 0.05, 0.05, 0.05)
    }else if( str_detect(mod_name, '2')){ 
      # HOI TYPE 2 
      init <- c(20, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30)
      lower <- c(1, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30)
      upper <- c(1e2, 4, 4, 4, 2, 2, 2)
    }else{
      # HOI TYPE 1 
      init <- c(20, -1, 0, 0, 0, 0, 0, 0)
      lower <- c(1, -2, 0, 0, 0, -1e-4, -1e-4, -1e-4)
      upper <- c(1e2, 0, 1e2, 1e2, 1e2, 1, 1, 1)
    }
  }else{
    if( str_detect(mod_name, '4')){ 
      # BASIC TYPE 4 
      init <- c(20, rep(1,3), rep(0, 3))
      lower <- c(1, rep(0,3), rep(0, 3))
      upper <- c(1e2, rep(2, 6))
    }else if( str_detect(mod_name, '3')){ 
      # BASIC  TYPE 3 
      init <- c(20, 1, 1, 1)
      lower <- c(1, 0, 0, 0)
      upper <- c(1e3, 4, 4, 4)
    }else if( str_detect(mod_name, '2') ){ 
      # BASIC  TYPE 2 
      init <- c(20, 1e-30, 1e-30, 1e-30)
      lower <- c(1, 1e-30, 1e-30, 1e-30)
      upper <- c(1e2, 4, 4, 4)
    }else{
      # BASIC  TYPE 1 
      init <- c(20, -1, 0, 0, 0)
      lower <- c(1, -2, 0, 0, 0)
      upper <- c(1e2, 0, 1e2, 1e2, 1e2)
    }
  }
  
  mm <- model.matrix(form, dat)
  
  y <- dat$fecundity
  init[1] <- max(y)
  
  fit <- optim(par = init, 
               fn = model, 
               mm = mm, 
               y = y, 
               method = 'L-BFGS-B', 
               sd = sd_grad[1], 
               lower = lower, 
               upper = upper) 
  
  for(j in 1:max_refit){ 
    
    fit <- optim(par = fit$par, 
                 fn = model, 
                 mm = mm, 
                 y = y, 
                 method = 'L-BFGS-B', 
                 sd = sd_grad[j], 
                 lower = lower, 
                 upper = upper) 
  }
  
  if(fit$convergence != 0){
    stop('fit1 did not converge')
  }
  
  return(fit)
}

prep_data <- function( sp, dat){ 
  focal <- paste0('F', sp)
  
  dat <- 
    dat %>% 
    filter_(.dots = paste0( 'focal == "', focal, '"')) %>% 
    distinct(id, focal, fecundity, comp_n, N1, N2, N3)
  
  return(dat)
}

predict_fit <- function(pars, dat, mod_name, form = form1){
  
  model <- eval(parse(text = mod_name))
  
  mm <- model.matrix(form, dat)
  y <- dat$fecundity
  
  y_pred <- model(pars, y = y, mm = mm, predict = T)
  
  return(y_pred)
}


fit_nls <- function(data, form, init_vals, max_comp = 2, ...){
  
  temp <- 
    data %>% 
    filter(comp_n < max_comp)  %>%
    mutate( y = fecundity) 
  
  test_fit <- nls(form, data = temp, start = init_vals, ... )
  
  return(test_fit ) 
  
} 