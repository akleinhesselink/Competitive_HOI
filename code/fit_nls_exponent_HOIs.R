rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')
library(mgcv)
nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'

dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter(comp_n < 3)

dat$y <- log(dat$fecundity)
gamfitBASIC <- gam(y ~ s(N1) + s(N2) + s(N3) + s(I(N1^2)) + s(I(N2^2)) + s(I(N3^2)), dat = dat)
gamfitHOI <- gam(y ~ s(N1) + s(N2) + s(N3) + s(I(N1*N2)) + s(I(N1*N3)) + s(I(N2*N3)), data = dat)
AIC(gamfitBASIC, gamfitHOI)

dat$pred <- predict(gamfitHOI)
plot(dat$y, dat$pred)

dat$y <- max(dat$fecundity)/dat$fecundity
gamfitBASIC <- gam(y ~ s(N1) + s(N2) + s(N3), dat = dat)
gamfitHOI <- gam(y ~ s(N1) + s(N2) + s(N3) + s(I(N1*N2)) + s(I(N1*N3)) + s(I(N2*N3)), data = dat)
AIC(gamfitBASIC, gamfitHOI)

dat$pred <- predict(gamfitBASIC)
plot(dat$y, dat$pred)

dat$pred <- predict(gamfitHOI)
plot(dat$y, dat$pred)



testfit <- 
  nls( y ~ beta + a[1]*N1 + a[2]*N2 + a[3]*N3, 
     data = dat, 
     start = list(beta = 0, a = c(0,0,0)))

plot( dat$y, predict(testfit))


testfit <- 
  nls( y ~ beta + a[1]*N1 + a[2]*N2 + a[3]*N3 + b[1]*I(N1*N2) + b[2]*I(N1*N3) + b[3]*I(N2*N3), 
       data = dat, 
       start = list(beta = 0, a = c(0,0,0), b = c(0,0,0)))

testfit
plot( dat$y, predict(testfit))




dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter(comp_n < 3)

testfit <- nls(log(fecundity) ~ log(lambda) - log( beta + 
                                                     N1^a[1] +  
                                                     N2^a[2] + 
                                                     N3^a[3] +
                                                     I(N2*N3)*b[1] +
                                                     I(N3*N1)*b[2] + 
                                                     I(N2*N1)*b[3]), 
               data = dat, 
               start = list(lambda = 44, beta = 1, a = c(0.2, 0.2, 0.9), b = c(0, 0, -0.001) ), 
               algorithm = 'port',
               lower = c(0, 0, 0, 0, 0, 0, 0, -1), 
               upper = c(500, 200, 1, 1, 1, 1, 1,  1))
testfit

dat$pred <- exp(predict(testfit))
dat %>% ggplot(aes( y = fecundity, x = N3)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N1,N2))), color = 'red')

dat %>% ggplot(aes( y = fecundity, x = N2)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N1,N3))), color = 'red')

dat %>% ggplot(aes( y = fecundity, x = N1)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = factor(paste(N2,N3))), color = 'red')



# try with species 3 left out 

dat <- prep_data(3, dat = readRDS(data_file))
dat <- 
  dat %>% 
  filter(comp_n < 3, N2 > 3 | N1 > 3 )

testfit <- nls(fecundity ~ lambda/( 1 + N1^(a[1]) + N2^(a[2]) ), 
             dat = dat, 
             start = list(lambda = 44, a = c(0.1, 0.1)), 
             algorithm = 'port', 
             lower = c(0,0,0), 
             upper = c(50, 5,5))

testfit <- nls(fecundity ~ lambda*((1 + N1^(a[1] + a[3]*N2) +  N2^(a[2] + a[4]*N1))^-1), 
               data = dat, 
               start = list(lambda = 44, a = c(0.1, 0.1, 0, 0)), 
               algorithm = 'port',
               lower = c(0, 0, 0, 0, 0), 
               upper = c(100, 1, 1,  1, 1))

testfit
coef(testfit)

dat$pred <- predict(testfit)

dat %>% 
  ggplot(aes( y = fecundity, x = N1)) + 
  geom_point() + 
  geom_line(aes(y = pred, group = N2), color = 'red')


