
rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

results_file <- 'data/ann_plant_sim1.rds'
pars_file    <- 'data/ann_plant_pars1.rds'

# set parameters ------------------------------------- 
nspp <- 3 
alphas <- matrix( c(1, 0.5, 0.1,  
                    0.3, 1, 0.4, 
                    0.1, 0.5, 1), nspp, nspp, byrow = T)

betas <- matrix(c(0.00,  0.01,  -0.01, 
                  0.00,  0.1,  0.00, 
                  0.05,  -0.01,  0.01), nspp, nspp, byrow = T)

lambdas <- c(24, 32, 41)
taus <- c(-1, -0.9, -0.7)

# 
maxdens <- 20
base <- 1.2

# -----------------------------------------------------
experiments <- make_experiments(maxdens, base, nspp)

out <- experiments
mm <- model.matrix(formHOI, experiments)

for( i in 1:nspp) { 
  out[,i] <- mod_bh_ll(pars = c(lambdas[i], taus[i], alphas[i, ], betas[i, ]), y = NA, mm = mm, predict = T)
}

names(out)[1:nspp] <- paste0('F', 1:nspp)

results <- 
  left_join(experiments, out, by = 'id') %>% 
  gather( focal, fecundity, starts_with('F')) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  group_by( id, focal ) %>% 
  mutate( comp_n = sum(density > 0)) %>% 
  filter( comp_n == 0 | density > 0 ) %>% 
  spread( competitor, density, fill = 0) %>% 
  ungroup()

results$focal_label <- paste0( 'focal\n', str_replace( results$focal, 'F', 'N'))

data.frame( species = paste0('N', 1:nspp), 
            lambda = lambdas, 
            alpha = alphas, 
            betas = betas, 
            tau = taus ) 
  
ann_plant_pars <- 
  data.frame( species = paste0('N', 1:nspp), 
              lambda = lambdas, 
              alpha = alphas, 
              betas = betas, 
              tau = taus ) %>% 
  gather( par, value, lambda:tau) %>%
  mutate( par = str_replace(par, '\\.', str_extract(species, '\\d+'))) %>%
  mutate( type = 'original')

saveRDS(ann_plant_pars, file = pars_file)
results_ann_plant <- results 
saveRDS(results_ann_plant, file = results_file)
