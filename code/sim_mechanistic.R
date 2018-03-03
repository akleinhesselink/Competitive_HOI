library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(stringr)
library(parallel)
rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
load('data/eqs.rda')
load('data/parms.rda')

# -------- simulate annual plant experiments -------------------------- # 

max_eqs <- apply( eqs, 2, max)
test <- eqs
test[test==0] <- NA 
min_eqs <- apply(test, 2, min, na.rm = T)

rs <- expand.grid( lapply( as.list(log(max_eqs*2, base = 2)), function(x) c(0, 2^(seq(0, x, 0.5))) ))

rs <- rs[-1, ]
rs <- rs %>% distinct()

rs <- split( rs, 1:nrow(rs))
rs_results <- mclapply( rs, run_experiment, parms, mc.cores = 4)

save(rs_results, file = 'data/rs_results.rda')
save(rs, file = 'data/rs.rda')
