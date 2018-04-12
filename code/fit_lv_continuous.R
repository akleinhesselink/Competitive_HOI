rm(list = ls())
library(deSolve)
source( 'code/LV_functions.R')
# set parameters ----------------------------------------- # 
nspp <- 2
alpha <- matrix(c(0.08, 0.001,
                  0.00, 0.09), nspp, nspp, byrow = T)

r <- c(0.2, 0.1)
parms <- list(r = r, alpha = alpha)
gen_time <- 2
times <- 500
# ------------------------------------------------------- # 
N_init <- c(1,1)
N_out <- matrix(NA, times, nspp )
N_out[1, ] <- N_init 
colnames(N_out) <- paste0('N', 1:nspp)
for( i in 2:times){ 
  N_out[i, ] <- run_lv_continuous(N = N_out[ i-1, ], model_fun = lv_continuous, parms, gen_time = gen_time)
}
matplot( N_out , type  = 'l') 

# run density gradient 
low <- 0 
high <- 20
by <- 2 

N <- as.matrix( expand.grid( rep(list(seq(low, high, by = by)), nspp) ))
colnames(N) <- paste0('N', 1:nspp)
N_out <- N
colnames(N_out) <- paste0('Y', 1:nspp)

for( i in 1:nrow(N)){ 
  N_out[i, ] <- run_lv_continuous(N = N[i, ], model_fun = lv_continuous, parms, gen_time = gen_time)
}

dat <- 
  data.frame( cbind(N, N_out/N)) %>% 
  gather( species, y, starts_with('Y')) %>% 
  filter( !is.na(y) ) 

temp <- split(dat, f = dat$species)

form <- make_form(nspp)
fit <- lapply( temp, fit_lv_discrete, form = form, par_init = c(0.3, 0.01, 0.01), lowers = c(0.01, -1, -1), uppers = c(2, 1, 1))

par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r
r_est
alpha_est <- par_ests[, -1]
alpha
round(alpha_est , 3)

# now try with HOI's  ---------------------------------- # 
form <- make_form(nspp, HOI = T)

fit <- lapply( temp, fit_lv_discrete, form = form, par_init = c(0.1, 0, 0, 0), lowers = c(0.01, -1, -1, -1), uppers = c(2, 1, 1, 1))

par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r
r_est
alpha_est <- par_ests[, -1]
alpha
round(alpha_est , 5)

# three species -------------------------------------------------- # 
nspp <- 3
r <- c(0.2, 0.1, 0.3)
alpha <- matrix( c(0.01,  0.02, 0.0, 
                   0.000, 0.01, 0.02, 
                   0.02,  0.00, 0.02), nspp, nspp, byrow = T)

parms <- list(r = r, alpha = alpha)
gen_time <- 2
times <- 500
# ------------------------------------------------------- # 
N_init <- c(2,1,3)
N_out <- matrix(NA, times, nspp )
N_out[1, ] <- N_init 
colnames(N_out) <- paste0('N', 1:nspp)
for( i in 2:times){ 
  N_out[i, ] <- run_lv_continuous(N = N_out[ i-1, ], model_fun = lv_continuous, parms, gen_time = gen_time)
}
matplot( N_out , type  = 'l') 

# run density gradient 
low <- 0 
high <- 50
by <- 2 

N <- as.matrix( expand.grid( rep(list(seq(low, high, by = by)), nspp) ))
colnames(N) <- paste0('N', 1:nspp)
N <- N[ apply( N, 1, function(x) sum(x > 0)) == 2, ]  # use only two species communities 

N_out <- N
colnames(N_out) <- paste0('Y', 1:nspp)

for( i in 1:nrow(N)){ 
  N_out[i, ] <- run_lv_continuous(N = N[i, ], model_fun = lv_continuous, parms, gen_time = gen_time)
}

dat <- 
  data.frame( cbind(N, N_out/N)) %>% 
  gather( species, y, starts_with('Y')) %>% 
  filter( !is.na(y) ) 

temp <- split(dat, f = dat$species)

form <- make_form(nspp)

fit <- lapply( temp, 
               fit_lv_discrete, 
               form = form, 
               par_init = c(0.3, 0.01, 0.01, 0.01),
               lowers = c(0.1, -0.5, -0.5, -0.5), 
               uppers = c(1, 0.5, 0.5, 0.5))

par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r_est
r
alpha_est <- par_ests[, -1]
alpha
round(alpha_est , 5)

# now try with HOIs ------------------------------------- # 
form <- make_form(nspp, HOI = T)
fit <- lapply( temp, 
               fit_lv_discrete, 
               form = form, 
               par_init = c(0.3, 0.01, 0.01, 0.01, 0, 0, 0), 
               lowers = c(0.05, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5), 
               uppers = c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))

par_ests <- do.call(rbind, lapply(fit, function(x) x$par))
r_est <- par_ests[, 1]
r
r_est
alpha_est <- par_ests[, -1]
alpha
round(alpha_est , 5)

form
