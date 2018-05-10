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

basic <- 'fecundity ~ beta0. + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3'
HOI <- 'fecundity ~ beta0. + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 +
                      gamma.[1]*I(N1^2) + gamma.[2]*I(N2^2) + gamma.[3]*I(N3^2) + 
                      beta.[1]*I(N1*N2) + beta.[2]*I(N1*N3) + beta.[3]*I(N2*N3)'

HOI2 <- 'fecundity ~ beta0. + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 +
                      gamma.[1]*I(N1^2) + gamma.[2]*I(N2^2) + gamma.[3]*I(N3^2) +
                      gamma.[4]*I(N1^3) + gamma.[5]*I(N2^3) + gamma.[6]*I(N3^3) + 
                      beta.[1]*I(N1*N2) + beta.[2]*I(N1*N3) + beta.[3]*I(N2*N3) + 
                      beta.[4]*I(N1^2*N2) + beta.[5]*I(N1^2*N3) + beta.[6]*I(N2^2*N3) + 
                      beta.[7]*I(N1*N2^2) + beta.[8]*I(N1*N3^2) + beta.[9]*I(N2*N3^2)'
  
lower_basic <- -100
upper_basic <- 1000
lower_HOI   <- -100
upper_HOI   <- 1000


fits <- plot <- list()

i <- 1
for ( i in 1:nspp){
  
  focal <- paste0('F', i)
  
  dat <- prep_data(i, dat = readRDS(data_file))
  
  dat <- dat %>% filter( comp_n  < 3)
  
  start_basic <- list(beta0. = 1, alpha. = rep(1,3))
  
  fit_basic <- nls(basic, 
                       data = dat, 
                       start = start_basic, 
                       algorithm = 'port', 
                       lower = lower_basic, 
                       upper = upper_basic)
  
  
  start_HOI <- list( beta0. = 1, alpha. = rep(1,3), gamma. = rep(1,6), beta. = rep(0, 9))
  
  fit_HOI   <- nls(HOI2, 
                   data = dat, 
                   start = start_HOI, 
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


