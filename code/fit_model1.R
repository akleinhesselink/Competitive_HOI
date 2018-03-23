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

for ( i in 1:nspp){
  
  focal <- paste0('F', i)
  
  dat <- prep_data(i, dat = readRDS(data_file))
  
  fit1 <- fit_model(dat = dat,
                    form = form1,
                    mod_name = model,
                    start_sd = start_sd,
                    min_sd = min_sd,
                    max_refit = max_refit)
  
  fitHOI <- fit_model(dat = dat,
                      form = formHOI,
                      mod_name = model,
                      start_sd = start_sd,
                      min_sd = min_sd,
                      max_refit = max_refit)
  
  dat$basic <- predict_fit(fit1$par, dat = dat, mod_name = model, form = form1)
  dat$HOI   <- predict_fit(fitHOI$par, dat = dat, mod_name = model, form = formHOI)
  
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



