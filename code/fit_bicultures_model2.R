rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3

init1 <- c(1, 1e-30, 1e-30, 1e-30)
initHOI <- c(1, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30)
lower1 <- c(1, 1e-30, 1e-30, 1e-30)
upper1 <- c(1e3, 1e2, 1e2, 1e2)
lowerHOI <- c(1, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30)
upperHOI <- c(1e3, 1e2, 1e2, 1e2, 1, 1, 1)

fits <- plot <- list()
for( i in 1:nspp){
  
  focal <- paste0('F', i)

  dat <- 
    readRDS('data/mechanistic_sim_bicultures.rds') %>% 
    filter_(.dots = paste0( 'focal == "', focal, '"')) %>% 
    distinct(id, focal, fecundity, N1, N2, N3)
  
  mm1 <- model.matrix(form1, dat)
  mmHOI <- model.matrix(formHOI, dat)
  
  y1 <- dat$fecundity
  
  fit1 <- optim(par = init1, 
                    fn = mod_bh2_ll, 
                    mm = mm1, 
                    y = y1, 
                    method = 'L-BFGS-B', 
                    sd = 0.5, 
                    lower = lower1, 
                    upper = upper1) 
  
  fitHOI <- optim(par = initHOI, 
                fn = mod_bh2_ll, 
                mm = mmHOI, 
                y = y1, 
                method = 'L-BFGS-B', 
                sd = 0.5, 
                lower = lowerHOI, 
                upper = upperHOI) 
  
  dat$basic <- mod_bh2_ll(fit1$par, y = y1, mm = mm1, predict = T)
  dat$HOI <- mod_bh2_ll(fitHOI$par, y = y1, mm = mmHOI, predict = T)
  
  dat <- dat %>% 
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
    theme(legend.position = c(0.85, 0.75), legend.box.just = 'right')
  
  plot[[i]] <- do.call( function(...) grid.arrange (..., ncol = nspp), list(p1, p2, p3) )
  fits[[i]] <- list(fit1 = fit1, fitHOI = fitHOI)
}

