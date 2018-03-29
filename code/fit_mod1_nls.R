rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

maxiter <- 10
count <- 0 
nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'
model <- "mod_bh_ll"

basic <- 'fecundity ~ lambda.*( (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau.) '
HOI <- 'fecundity ~ lambda.*( (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 + betas.[1]*I(N1*N2) + betas.[2]*I(N1*N3) + betas.[3]*I(N2*N3))^tau.) '

lower_basic <- c(1, -1.5, rep(0, 3))
upper_basic <- c(1e3, -0.1, rep(100, 3))
lower_HOI   <- c(lower_basic, rep(-0.2, 3))
upper_HOI   <- c(upper_basic, rep(50, 3))


fits <- plot <- list()

i <- 3
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

coef(fits[[3]]$basic)
coef(fits[[3]]$HOI)

sum(resid(fits[[3]]$basic)^2)
sum(resid(fits[[3]]$HOI)^2)

HOI <-  data.frame( species = c('N1', 'N2', 'N3'), 
                    do.call(rbind, lapply( fits, function(x) coef(x$HOI))))

fitted <- HOI %>% 
  gather( par, value, - species ) %>% 
  mutate( par  = str_replace(par, '\\.', (str_extract(species, '\\d+')))) %>% 
  mutate( type = 'fitted')

compare_parameters('data/ann_plant_pars1.rds', fitted )


dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter( N1 == 0 , N2 ==0 ) 
testfit3 <- nls(fecundity ~ lambda*((1 + N3^a[1])^tau), data = dat, start = list(lambda = 10, a = c(1), tau = -1)  )
dat$pred <- predict(testfit3)
dat %>% ggplot(aes( y = fecundity, x = N3)) + 
  geom_point() + 
  geom_line(aes(y = pred), color = 'red')

dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter( N3 == 0, N2 ==0 ) 
testfit1 <- nls(fecundity ~ lambda*((1 + N1^a[1])^tau), data = dat, start = list(lambda = 40, a = c(1), tau = -1)  )
dat$pred <- predict(testfit1)
dat %>% ggplot(aes( y = fecundity, x = N1)) + 
  geom_point() + 
  geom_line(aes(y = pred), color = 'red')

dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter( N3 == 0, N1 ==0 ) 
testfit2 <- nls(fecundity ~ lambda*((1 + N2^a[1])^tau), data = dat, start = list(lambda = 40, a = c(1), tau = -1)  )
dat$pred <- predict(testfit2)
dat %>% ggplot(aes( y = fecundity, x = N2)) + 
  geom_point() + 
  geom_line(aes(y = pred), color = 'red')



dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter(comp_n < 3, N2 == 0 )

testfit <- nls(fecundity ~ lambda*((1 + N1^(a[1] + a[3]*N3) +  N3^(a[2]))^-1), 
               data = dat, 
               start = list(lambda = 44, a = c(0.2, 0.9, 0)), 
               algorithm = 'port',
               lower = c(0, 0, 0, -30), 
               upper = c(50, 5, 5,  4))

coef(testfit)
dat$pred <- predict(testfit)
dat %>% ggplot(aes( y = fecundity, x = N3)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N1,N2))), color = 'red')


dat %>% ggplot(aes( y = fecundity, x = N2)) + 
  geom_point() + 
  geom_line(aes(group = factor(paste(N1,N3))), color = 'blue', alpha = 0.1) + 
  geom_line(aes(y = pred, group = factor(paste(N1,N3))), color = 'red')


#
dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter(comp_n == 2)

testfit <- nls(fecundity ~ lambda*((1 + N1^a[1] + N2^a[2] +  N3^a[3])^tau), 
               data = dat, 
               start = list(lambda = 44, a = c(0.2, 0.2, 0.9), tau = -0.5), 
               algorithm = 'port',
               lower = c(0, 0, 0, 0, -1.1), 
               upper = c(50, 5, 5, 5,  -0.1))

coef(testfit)
dat$pred <- predict(testfit)
dat %>% ggplot(aes( y = fecundity, x = N3)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N1,N2))), color = 'red')

dat %>% ggplot(aes( y = fecundity, x = N2)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N1,N3))), color = 'red')

dat %>% ggplot(aes( y = fecundity, x = N1)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N2,N3))), color = 'red')

dat1 <- dat %>% filter( N3 > 0 & N2 > 0)

testfit <- nls(fecundity ~ 48.5*((1 + 
                                      N1^0.12 + 
                                      N2^0.3 +  
                                      N3^0.99 + 
                                      I(N3*N2)^a[1])^1.1), 
               data = dat1, 
               start = list(a = 0.2), 
               algorithm = 'port',
               lower = c(-30), 
               upper = c(0.3))

coef(testfit)

dat1 <- dat1 %>% mutate( I1 =  N1*N2)
dat1$I1
