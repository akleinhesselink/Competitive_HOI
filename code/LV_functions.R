
lv_continuous <- function(t, N, parms){ 
  
  with(parms, { 
    
    dN <-  r*N*( 1 - alpha %*% N )
    
    dN <- unlist(dN)
    
    return(list(dN))
  })
  
}

run_lv_continuous <- function( N, model_fun, parms, gen_time ) { 
  
  out <- ode(N, seq(1, gen_time, by = 0.1), model_fun, parms = parms)
  
  tail( out, 1)[-1]
  
}



mod_lv_discrete <- function(par, y, mm, sd = 1, predict = FALSE){ 
  
  r <- par[1]
  alpha <- par[-1]
  
  mu <- 1 + r*(1 - mm %*% alpha)
  
  neg_ll <- sum( - dnorm( mu, y, sd, log = T) )
  
  if(predict){ 
    return(mu)
  }else if(!predict){ 
    return(neg_ll) 
  }
}

make_form <- function(nspp, HOI = F){ 
  
  spp <- paste0( 'N', 1:nspp)
  
  form <- paste0(spp, collapse = ' + ')
  
  if( HOI ) { 
    hois <- apply( t( combn(spp, 2)), 1, paste0, collapse = '*')
    form <- paste(form, paste0( 'I(', paste0(hois, collapse = ') + I('), ')'), sep = ' + ')
  }
  
  form <- as.formula( paste0(' ~ -1 + ', form) )
  
  return(form)
}

run_time_series <- function(nspp, form, N_init, times, r, alpha){ 
  
  N_out <- data.frame( matrix(NA, times, nspp ) )
  N_out[1, ] <- N_init 
  colnames(N_out) <- paste0('N', 1:nspp)
  
  for( i in 2:times){ 
    for(j in 1:nspp){ 
      mm <- model.matrix(form, N_out[i-1, ])
      N_out[i, j] <- N_out[i-1, j]*mod_lv_discrete(par = c(r[j], alpha[j,]), y = NA, mm = mm, predict = T)  
    }
  }
  return(N_out)
}

run_density_gradient <- function(nspp, form, low, high, by, r, alpha){ 
  
  N <- data.frame( as.matrix( expand.grid( rep(list(seq(low, high, by = by)), nspp) ) ) )
  colnames(N) <- paste0('N', 1:nspp)
  N_out <- N
  colnames(N_out) <- paste0('Y', 1:nspp)
  
  for( i in 1:nrow(N)){ 
    for(j in 1:nspp){ 
      mm <- model.matrix(form, N[i,])
      N_out[i, j] <- mod_lv_discrete(par = c(r[j], alpha[j,]), y = NA, mm = mm, predict = T)  
    }
  }
  dat <- 
    cbind( N, N_out) %>% 
    gather( species, y, starts_with('Y')) %>% 
    filter( !is.na(y) ) 
  
  return(dat)
}

fit_lv_discrete <- function(form, temp, par_init, lowers, uppers){
  
  mm <- model.matrix(form, temp) 
  fit <- optim(par = par_init, 
               mod_lv_discrete, 
               y = temp$y, 
               mm = mm, 
               method = 'L-BFGS-B', 
               lower = lowers, 
               upper = uppers)
  return(fit)
}




