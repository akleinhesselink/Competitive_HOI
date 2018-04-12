rm(list = ls())
library(deSolve)
library(tidyverse)
source('code/LV_functions.R')

# set parameters ------------------------------------------------- # 
nspp <- 2
alpha <- matrix(c(0.02, 0.01, 
                  0.00, 0.08), nspp, nspp, byrow = T)

r <- c(0.15, 0.2)
parms <- list(r = r, alpha = alpha)

# run two species simulation ------------------------------------------------- # 
nspp <- 2
N_init <- rep(1, nspp)
times <- 100

form <- make_form(nspp)
N_out <- run_time_series(nspp = 2, form = form, N_init = c(1,1), times = 100, r = r, alpha = alpha)

matplot( N_out , type = 'l')

# density gradient experiments --------------------------------------------- # 
dat <- run_density_gradient(2, form, low = 0, high = 20, by = 1, r = r, alpha = alpha )
temp <- split(dat, f = dat$species)

fit <- lapply( temp, 
               fit_lv_discrete, 
               form = form, 
               par_init = c(0.1, 0, 0), 
               lowers = c(0.01, 0, 0), 
               uppers = c(1, 1, 1))

par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r_est
r
alpha_est <- par_ests[, -1]
alpha == round( alpha_est, 3)
alpha 
round( alpha_est, 3)

# Three species experiment -------------------------------------------------- # 
nspp <- 3
r <- c(0.1, 0.2, 0.3)
alpha <- matrix( c(0.02, 0.01, 0.005, 
                   0.005, 0.03, 0.01, 
                   0.008, 0.01, 0.02), nspp, nspp, byrow = T)

form <- make_form(nspp )
out <- run_time_series(nspp, form, N_init = c(1,1,1), times = 100, r = r, alpha = alpha )
matplot(out, type = 'l')

# density gradient
dat <- run_density_gradient(nspp, form, low = 0, high = 20, by = 2, r = r, alpha = alpha)

temp <- split(dat, f = dat$species)

fit <- lapply( temp, 
               fit_lv_discrete, 
               form = form, 
               par_init = c(0.1, 0, 0, 0), 
               lowers = c(0.01, 0, 0, -0.5), 
               uppers = c(1, 1, 1, 0.5))


par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r
r_est
alpha_est <- par_ests[, -1]
alpha
round(alpha_est , 3)
alpha == round( alpha_est, 3)

# three species with HOIs ------------------------------------------------- # 
nspp <- 3
r <- c(0.1, 0.2, 0.3)
alpha <- matrix( c(0.02, 0.01, 0.005, 
                   0.005, 0.03, 0.01, 
                   0.008, 0.01, 0.02), nspp, nspp, byrow = T)

# HOI matrix 
beta <- matrix( c(0, 0, -0.001,  
                  0, 0, 0, 
                  0.001, 0, 0), nspp, nspp, byrow = T ) 

alpha <- cbind(alpha, beta)                

form <- make_form(nspp, HOI = T )
out <- run_time_series(nspp, form, N_init = c(1,1,1), times = 100, r = r, alpha = alpha )
matplot(out, type = 'l')

# density gradient 
dat <- run_density_gradient(nspp, form, low = 0, high = 40, by = 2, r = r, alpha = alpha)

temp <- split(dat, f = dat$species)

fit <- lapply( temp, 
               fit_lv_discrete, 
               form = form, 
               par_init = c(0.3, 0, 0, 0, 0, 0, 0), 
               lowers = c(0.01, 0, 0, 0, -0.5, -0.5, -0.5), 
               uppers = c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))


par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r_est
r
alpha_est <- par_ests[, -1]
alpha
round(alpha_est , 3)
alpha == round( alpha_est, 4)


