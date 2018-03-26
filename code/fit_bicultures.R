rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3

init1 <- c(1, -1, 0, 0, 0)
initHOI <- c(1, -1, 0, 0, 0, 0, 0, 0)
lower1 <- c(1, -10, 0, 0, 0)
upper1 <- c(1e3, 0, 1e2, 1e2, 1e2)
lowerHOI <- c(1, -10, 0, 0, 0, 0, 0, 0)
upperHOI <- c(1e3, 0, 1e2, 1e2, 1e2, 1, 1, 1)

focal <- 'F3'

dat <- 
  readRDS('data/mechanistic_sim_bicultures.rds') %>% 
  filter_(.dots = paste0( 'focal == "', focal, '"')) %>% 
  distinct(id, focal, fecundity, comp_n, N1, N2, N3)

mm1 <- model.matrix(form1, dat)
mmHOI <- model.matrix(formHOI, dat)

y1 <- dat$fecundity

fit1 <- optim(par = init1, 
                  fn = mod_bh_ll, 
                  mm = mm1, 
                  y = y1, 
                  method = 'L-BFGS-B', 
                  sd = 1, 
                  lower = lower1, 
                  upper = upper1) 

fitHOI <- optim(par = initHOI, 
              fn = mod_bh_ll, 
              mm = mmHOI, 
              y = y1, 
              method = 'L-BFGS-B', 
              sd = 1, 
              lower = lowerHOI, 
              upper = upperHOI) 

dat$basic <- mod_bh_ll(fit1$par, y = y1, mm = mm1, predict = T)
dat$HOI <- mod_bh_ll(fitHOI$par, y = y1, mm = mmHOI, predict = T)

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

do.call( function(...) grid.arrange (..., ncol = nspp), list(p1, p2, p3) )














