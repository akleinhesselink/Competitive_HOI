rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'
model <- "mod_bh_ll"

fits <- plot <- list()
i <- 3
for ( i in 1:nspp){
  
  focal <- paste0('F', i)
  
  dat <- prep_data(i, dat = readRDS(data_file))
  
  dat <- do.call(rbind, list(dat,dat,dat,dat,dat))
  dat$y <- dat$fecundity
  
  fit1 <- nls(y ~ l*((1 + a1*N1 + a2*N2 + a3*N3)^ta), data = dat, start = c(l = max(dat$y), a1 = 1, a2 = 1, a3 = 1, ta = -1))
  
  fit1 <- nls(y ~ l*(1 + a1*N1 + a2*N2 + a3*N3)^ta, data = dat , start = coef(fit1))
  dat$basic <- predict(fit1)
  
  fitHOI <- nls(y ~ l*(1 + a1*N1 + a2*N2 + a3*N3 + b1*sqrt(I(N1*N2)) + b2*sqrt(I(N1*N3)) + b3*sqrt(I(N3*N2)))^ta, 
                       data = dat, 
                       start = c(coef(fit1), b1 = 1, b2 = 1, b3 = 1))
  
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
