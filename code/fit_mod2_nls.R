rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'

basic <- 'fecundity ~ lambda./(1 + N1^alpha.[1] + N2^alpha.[2] + N3^alpha.[3])'
HOI <- 'fecundity ~ lambda./(1 + N1^(alpha.[1] + betas.[1]*N2 + betas.[2]*N3)  + N2^(alpha.[2] + betas.[3]*N1 + betas.[4]*N3) + N3^(alpha.[3] + betas.[5]*N1 + betas.[6]*N2) )'

curve(3^(1/(1 + x)), -2, 10)

lower_basic <- c(1, rep(0, 3))
upper_basic <- c(1e3, rep(2, 3))
lower_HOI   <- c(lower_basic, rep(0, 6))
upper_HOI   <- c(upper_basic, rep(1, 6))


fits <- plot <- list()

i <- 3
for ( i in 1:nspp){
  
  focal <- paste0('F', i)
  
  dat <- prep_data(i, dat = readRDS(data_file))
  
  dat <- dat %>% filter(N1 > 0 | N2 > 0 | N3 > 0)
  
  start_basic <- list(lambda. = max(dat$fecundity), alpha. = rep(0.5,3))
  
  fit_basic <- try(nls(basic, 
                       data = dat, 
                       start = start_basic, 
                       algorithm = 'port', 
                       lower = lower_basic, 
                       upper = upper_basic), silent = T)
  
  
  while( as.character(class(fit_basic)) == 'try-error'){ 
    
    if(start_basic$alpha.[1] > upper_basic[2]){ stop('tau too high')}
    
    start_basic$alpha.[1] <- start_basic$alpha.[1] + 0.01
    
    fit_basic <- 
      try(nls(basic, 
              data = dat, 
              start = start_basic, 
              algorithm = 'port', 
              lower = lower_basic, 
              upper = upper_basic), silent = T)
    
  }
  start_HOI <- list( lambda. = coef(fit_basic)[1], 
                     alpha. = coef(fit_basic)[2:4], 
                     betas. = c(0,0,0,0,0,0))

  fit_HOI   <- nls(HOI, 
                   data = dat, 
                   start = start_HOI, 
                   algorithm = 'port', 
                   lower = lower_HOI, 
                   upper = upper_HOI)
  fit_basic
  fit_HOI
  
  dat$basic <- predict(fit_basic)
  dat$HOI   <- predict(fit_HOI)
  
  dat <-
    dat %>%
    gather( type, pred, c(basic, HOI))
  
  p1 <- plot_two_sp(dat, focal = focal, C1 = 'N1', C2 = 'N2', C3 = 'N3') +
    geom_line(aes(y = pred, linetype = type)) +
    scale_linetype_manual(values = c(2,3), guide = F) +
    theme(legend.position = c(0.85, 0.8))
  
  p2 <- plot_two_sp(dat, focal = focal, C1 = 'N2', C2 = 'N3', C3 = 'N1') +
    geom_line(aes(y = pred, linetype = type)) +
    scale_linetype_manual(values = c(2,3), guide = F) +
    theme(legend.position = c(0.85, 0.8))
  
  p3 <- plot_two_sp(dat, focal = focal, C1 = 'N3', C2 = 'N1', C3 = 'N2') +
    geom_line(aes(y = pred, linetype = type)) +
    scale_linetype_manual(values = c(2,3)) +
    theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')
  
  plot[[i]] <- do.call( function(...) grid.arrange (..., ncol = nspp), list(p1, p2, p3) )
  fits[[i]] <- list(basic = fit_basic, HOI = fit_HOI)
  
}

