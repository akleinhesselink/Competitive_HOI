rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

maxiter <- 10
count <- 0 
nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/ann_plant_sim1.rds'
model <- "mod_bh_ll"

basic <- 'fecundity ~ lambda.*( (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau.) '
HOI <- 'fecundity ~ lambda.*( (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 + betas.[1]*I(N1*N2) + betas.[2]*I(N1*N3) + betas.[3]*I(N2*N3))^tau.) '

lower_basic <- c(1, -1.5, rep(0, 3))
upper_basic <- c(1e3, -0.5, rep(2, 3))
lower_HOI   <- c(lower_basic, rep(-0.1, 3))
upper_HOI   <- c(upper_basic, rep(0.1, 3))


fits <- plot <- list()


for ( i in 1:nspp){

  focal <- paste0('F', i)

  dat <- prep_data(i, dat = readRDS(data_file))
  
  start_basic <- list(lambda. = max(dat$fecundity), tau. = -1.5, alpha. = rep(1,3))
  
  fit_basic <- try(nls(basic, 
                 data = dat, 
                 start = start_basic, 
                 algorithm = 'port', 
                 lower = lower_basic, 
                 upper = upper_basic), silent = T)
  
  
  while( as.character(class(fit_basic)) == 'try-error'){ 
    
    if(start_basic$tau. > upper_basic[2]){ stop('tau too high')}
    
    start_basic$tau. <- start_basic$tau. + 0.01

    fit_basic <- 
      try(nls(basic, 
            data = dat, 
            start = start_basic, 
            algorithm = 'port', 
            lower = lower_basic, 
            upper = upper_basic), silent = T)

  }
  
  fit_HOI   <- nls(HOI, 
                 data = dat, 
                 start = list( lambda. = coef(fit_basic)[1], tau. = coef(fit_basic)[2], alpha. = coef(fit_basic)[3:5], betas. = rep(0,3)), 
                 algorithm = 'port', 
                 lower = lower_HOI)
  
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


HOI <-  data.frame( species = c('N1', 'N2', 'N3'), 
                    do.call(rbind, lapply( fits, function(x) coef(x$HOI))))

fitted <- HOI %>% 
  gather( par, value, - species ) %>% 
  mutate( par  = str_replace(par, '\\.', (str_extract(species, '\\d+')))) %>% 
  mutate( type = 'fitted')

compare_parameters('data/ann_plant_pars1.rds', fitted )


