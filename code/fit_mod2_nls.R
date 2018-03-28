rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'

fits <- plot <- list()
i <- 1
for ( i in 1:nspp){
  
  focal <- paste0('F', i)
  
  dat <- prep_data(i, dat = readRDS(data_file))
  
  dat <- do.call(rbind, list(dat,dat,dat,dat,dat))
  dat$y <- dat$fecundity
  
  nls_form1 <- 'y ~ l/(1 + N1^a1 + N2^a2 + N3^a3)'
  nls_formHOI <- 'y ~ l/(1 + N1^a1 + N2^a2 + N3^a3 + I(N1*N2)^b1 + I(N1*N3)^b2 + I(N2*N3)^b3)'
         
  fit1 <- nls(nls_form1, 
              data = dat, 
              start = c(l = max(dat$y), a1 = 1e-3, a2 = 1e-3, a3 = 1e-3), 
              algorithm = 'port', 
              lower = 1e-3, 
              upper = c(1e2, 2, 2, 2))
  
  fit1 <- nls(nls_form1, 
              data = dat, 
              start = coef(fit1), 
              algorithm = 'port',
              lower = 1e-30, 
              upper = c(1e2, 2,2,2))
  
  dat$basic <- predict(fit1)
  
  
  fitHOI <- nls( nls_formHOI, 
                data = dat, 
                start = c(coef(fit1), b1 = 0.01, b2 = 0.01, b3 = 0.01), 
                lower = 1e-10, 
                upper = c(1e2, rep(1, 6)), 
                algorithm = 'port', 
                control  = list(warnOnly = T) )
  
  coef(fitHOI)
  
  nls.control()
  
  fitHOI <- nls(y ~ l*(1 + a1*N1 + a2*N2 + a3*N3 + b1*sqrt(I(N1*N2)) + b2*sqrt(I(N1*N3)) + b3*sqrt(I(N3*N2)))^ta, 
                data = dat, 
                start = coef(fitHOI))
  
  dat$HOI <- predict(fitHOI)
  
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
  fits[[i]] <- list(fit1 = fit1, fitHOI = fitHOI)
  
}
