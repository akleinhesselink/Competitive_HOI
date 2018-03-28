rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/ann_plant_sim1.rds'
model <- "mod_bh_ll"


basic <- 'fecundity ~ lmbda/(1 + a[1]*N1 + a[2]*N2 + a[3]*N3)'
HOI <- 'fecundity ~ lmbda/(1 + a[1]*N1 + a[2]*N2 + a[3]*N3 + b[1]*I(N1*N2) + b[2]*I(N1*N3) + b[3]*I(N2*N3))'
lower_basic <- c(1, 0, 0, 0)
lower_HOI   <- c(lower_basic, -0.5, -0.5, -0.5)


fits <- plot <- list()

i <- 1
for ( i in 1:nspp){

  focal <- paste0('F', i)

  dat <- prep_data(i, dat = readRDS(data_file))
  
  start_basic <- list(lmbda = max(dat$fecundity), a = rep(0,3))
  
  fit_basic <- nls(basic, 
                 data = dat, 
                 start = start_basic, 
                 algorithm = 'port', 
                 lower = lower_basic)
  
  fit_HOI   <- nls(HOI, 
                 data = dat, 
                 start = list( lmbda = coef(fit_basic)[1], a = coef(fit_basic)[2:4], b = rep(0,3)), 
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
