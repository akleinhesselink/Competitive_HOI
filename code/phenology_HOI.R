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
times <- 125             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.9, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(98, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedlings <- c(1, 0, 0)      # number of seedlings 
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q)
rm( list = names(parms) )

# -------- simulate annual plant experiments -------------------------- # 
comp_grad <- c(0:90)  # number of competitors
experiments <- expand.grid(as.list(rep(list ( comp_grad), 3))) # response surface experiment 
experiments <- experiments[-1, ]
experiments <- experiments %>% distinct()
experiments <- experiments[ rowSums(experiments) < 30, ] 

fecundity <- phenology <- data.frame( matrix( NA, nrow = nrow(experiments), ncol = 3))
per_capita_use <- use <- out <- list()

run_experiment <- function(seedlings, soil_m, seedling_mass) { 
  seedlings <- as.numeric(seedlings)
  State <- c(soil_m, seedlings*seedling_mass)
  ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
}

start1 <- proc.time()
out <- mclapply( 
  split(experiments, f = 1:nrow(experiments)), 
  run_experiment, 
  soil_m, 
  seedling_mass, mc.cores = 4
  )
t1 <- proc.time() - start1
t1
use <- out

get_resource_use <- function(State, f, parms){ 
  State[3:5]*f(State[2], parms$r, parms$K)
}

use <- lapply( out, function(x) t(apply( x, 1, get_resource_use, f, parms)))
experiments <- split(experiments,1:nrow(experiments))
get_per_capita_use <- function(x, y) sweep( x, 2, as.numeric(y), '/')
per_capita_use <- mapply(use, experiments, FUN = get_per_capita_use, SIMPLIFY = F)
phenology <- lapply( out, function( x ) apply(x[, c(3:5)], 2, find_phenology))

##

experiments <- do.call(rbind, experiments)
max_biomass <- do.call( rbind, lapply( out, function( x ) apply( x[, c(3:5)], 2, max)))
fecundity <- (max_biomass*conversion)/seedling_mass
phenology <- do.call(rbind, phenology )

lambda <- diag( as.matrix( fecundity [ apply(experiments, 1, sum) == 1, ]  )) # calculate lambdas 
y <- fecundity/experiments                                                    # calculate fecundity in all experiments 
data <- data.frame(experiments, y, phenology/10)
names(data ) <- c('N1', 'N2', 'N3', 'Y1', 'Y2', 'Y3', paste0('PH', c(1:3)))

form1 <- paste('~ -1 + N1 + N2 + N3')
form2 <- paste(form1, '+ I(N1*N2) + I(N1*N3) + I(N2*N3)')

all_forms <- c(form1, form2)

gg_intra <- gg_inter <- gg_phenology <- list()
results <-  list()
i = 1

# fit annual plant model parameters --------------------------------------- # 
for(i in 1:3){ 
  # loop through each species and fit annual plant model parameters ------------------------------------------- # 
  data1 <- data
  data1[, i] <- data1[, i] - 1            # remove focal from competive neighborhood 
  data1$y  <- data1[, i + 3]              # focal per capita seed production
  data1$data <- as.matrix(data1[,c(1:3)])
  data2 <- data1[is.finite(data1$y),]  

  terms <- lapply( all_forms, function(x) length(str_split(x, pattern = '\\+')[[1]] ))
  forms <- lapply(all_forms, as.formula)
  
  pars <- lapply( terms, function(x, ...) {c( ..., rep(0, x - 1), -1 )} , Z = lambda[i] )  
  out <- mapply(par = pars, form = forms, FUN = optim, MoreArgs = list(data = data2, fn = mod_bh, method = 'BFGS', control = list( maxit = 1e9)), SIMPLIFY = F)
  results[[i]] <- out 
  
  pars <- lapply( out , function(x) x$par)
  predictions <- mapply(par = pars, form = forms, FUN = mod_bh2, MoreArgs = list(data = data1, predict = T))
  predictions <- data.frame(predictions)  
  names(predictions) <- paste0(paste0('pY', i, '.model'), 1:length(forms))
  data <- cbind(data, predictions)
} 

data <- 
  data %>% 
    gather(prediction, value,  starts_with('p', ignore.case = F)) %>% 
    gather( species, observed, Y1:Y3) %>% 
    separate( prediction, c('pred_species', 'model'), sep = '\\.') %>% 
    filter(str_extract(pred_species, 'Y\\d+') == str_extract(species, 'Y\\d+')) %>% 
    select( - pred_species) %>% 
    spread( model, value ) %>% 
    filter( is.finite(observed)) %>% 
  mutate( diff1 = (model1 - observed)/observed, diff2 = (model2 - observed)/observed) %>% 
  gather( model, deviation, diff1, diff2)

ggplot( data %>% filter( N1 ==  4), aes( x = N2, y = N3, fill = deviation )) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  facet_wrap(~ model)

ggplot( data %>% filter( N2 == 4), aes( x = N1, y = N3, fill = deviation )) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  facet_wrap(~ model)

ggplot( data %>% filter( N3 == 4), aes( x = N1, y = N2, fill = deviation )) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  facet_wrap(~ model)


basic_plot <- function(data, focal, hold_at = 4) { 
  N_focal <- paste0('N', focal)
  Y_focal <- paste0('Y', focal)
  comp <- paste0('N', c(1:3)[-focal])
  
  data <- 
    data %>% 
    filter( species == Y_focal) %>% 
    filter_(.dots = paste0(N_focal, ' == ', hold_at)) %>% 
    filter_(.dots = paste0(comp[1], ' %in% c(0, 1, 4, 8, 12)')) %>% 
    filter_(.dots = paste0(comp[2], ' %in% c(0, 1, 4, 8, 12)')) %>% 
    select( N1:N3, species, observed, model1, model2) %>% 
    gather(type, predicted, model1:model2)

  data[, comp[2]] <- as.factor(data[, comp[2]])
  
  ymx <- max( c(data$observed, data$predicted))
  
  data %>% 
    ggplot(aes_string( x = comp[1], y = 'observed', color = comp[2])) +
      geom_point() +
      geom_line( aes(y = predicted, linetype = type)) + 
      facet_wrap( ~ type ) + 
    ylim ( 0, ymx + 1)
}


data %>% filter( species == 'Y3', N3 == 4, N1 == 1, N2 == 1) %>% head()
data %>% basic_plot(focal = 3, hold_at = 9)

results

data1 <- 
  data1 %>% 
  select ( -data) %>% 
  gather( type, val, c(y, starts_with('type'))) 

data1$type <- factor(data1$type, labels = c('pred1', 'HOI', 'observed'))
data1$type_form <- factor(data1$type, labels = c('inter + intra', 'HOI', 'observed'))
data_intra <- data1[ rowSums( data1[, c(1:3)[-i]] ) == 0, ]

gg_intra[[i]] <- 
  ggplot( data_intra %>% filter( type != 'observed'), aes_string( x = names(data1)[i], y = names(data1)[i + 3])) + 
  geom_point( color = my_colors[i]) + 
  geom_line(aes( y = val), color = my_colors[i], linetype = 2) + 
  xlab( paste0( 'Number of intraspecific competitors')) +
  ylab( paste0( 'Per capita seed production')) +
  ggtitle( paste0('Species ', i)) + 
  my_theme + facet_wrap(~type_form)

data_inter <- data1[ data1[ , i ] == 5, ] 

data_inter$c1 <- data_inter[, c(1:3)[-i][1]]
data_inter$c2 <- as.factor(data_inter[, c(1:3)[-i][2]])

gg_inter[[i]] <- 
  ggplot(data_inter %>% filter( type != 'observed'), aes_string( x = 'c1', y = names(data_inter)[4:6][i], color = 'c2')) + 
  geom_point() + 
  geom_line(aes(y = val), linetype = 2) + 
  scale_color_discrete(paste( 'Density of', names(data_inter)[1:3][-i][2])) + 
  xlab( paste( 'Density of', names(data_inter)[1:3][-i][1])) + 
  ylab( paste( 'Per capita seed production of', names(data_inter)[i])) + 
  my_theme + 
  facet_wrap(~type_form)

pheno_data <- data1
sel_ph <- paste0( 'PH', i)
sel_self <- paste0( 'N', i)
sel_comp <- paste0( 'N', c(1:3)[-i])

intra_pheno <- pheno_data[ (pheno_data[,i] > 0 | rowSums(pheno_data[, c(1:3)]) == 0 ), ] %>% select_(sel_self, sel_ph) 
intra_pheno[, 1 ] <- intra_pheno[, 1] + 1 
intra_pheno <- intra_pheno %>% gather( competitor, density, 1)

inter_pheno <- pheno_data[ (pheno_data[,i] == 0 & rowSums(pheno_data[, c(1:3)[-i]]) > 0 ), ] %>% select_(sel_comp[1], sel_comp[2], sel_ph) 
inter_pheno <- inter_pheno %>% gather( competitor, density , 1:2 ) %>% filter( density > 0 )

pheno_data <- rbind( intra_pheno, inter_pheno)  
names(pheno_data)[1] <- 'day'

gg_phenology[[i]] <- 
  ggplot(pheno_data, aes( x = density, y = day, color = competitor )) + 
  geom_point() + 
  geom_line() + 
  xlab( 'Density of competitors') + 
  ylab( paste0( 'Day of flowering of N', i)) + 
  scale_color_discrete('Competitor Species') + 
  my_theme 



gg_intra[[1]] %+% subset( gg_intra[[1]]$data)
gg_intra[[2]] %+% subset( gg_intra[[2]]$data)
gg_intra[[3]] %+% subset( gg_intra[[3]]$data) 

# Compare outcomes in 2 and 3 species communities -------------------------------- # 
gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,6))
gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,6))
gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,6))


# compare error of fits 
fits <- data.frame( do.call( rbind, ( lapply( results, function(x) unlist(lapply(x, function(x) x$value))))))
fits$species <- paste('species', 1:3)
fits <- fits %>% gather( model, 'MSE', X1:X4, -species)
fits$model_label <- factor(fits$model, labels = levels( gg_inter[[1]]$data$type_form)[1:4])
fits <- fits %>% group_by( species ) %>% mutate( MSE_rel = MSE/max(MSE))

ggplot( fits %>% filter( model == 'X4' ), aes( x = species, y = MSE_rel)) + 
  geom_bar(stat = 'identity')

fit_plot <- ggplot( fits, aes( x = model_label, y = MSE) ) + 
  geom_bar(stat= 'identity') + 
  facet_grid(species ~. , scales = 'free') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5 ), axis.title.x = element_blank()) 

fit_plot

ggsave('figures/fit_plot.png', fit_plot, width = 5, height = 3.7)
ggsave('figures/intra1.png', gg_intra[[1]] %+% subset( gg_intra[[1]]$data, type %in% c('pred1', 'pred2')), width = 5, height = 3.3 ) 
ggsave('figures/intra2.png', gg_intra[[2]] %+% subset( gg_intra[[2]]$data, type %in% c('pred1', 'pred2')), width = 5, height = 3.3 ) 
ggsave('figures/intra3.png', gg_intra[[3]] %+% subset( gg_intra[[3]]$data, type %in% c('pred1', 'pred2')), width = 5, height = 3.3 ) 

ggsave( 'figures/inter1.png', gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,6)), width = 5.9, height = 3.3)
ggsave( 'figures/inter2.png', gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,6)), width = 5.9, height = 3.3)
ggsave( 'figures/inter3.png', gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,6)), width = 5.9, height = 3.3)

#
test <- results
names(test) <-  c('N1', 'N2', 'N3')
test <- lapply( test, function( x) {names(x) <- c('m1', 'm2'); x} )
result_dat <- data.frame( value = unlist( test), label = names(unlist(test)))

test <- 
  result_dat %>% 
  separate( label, c('species', 'model', 'par'), sep = '\\.', extra = 'warn', fill = 'left', convert = T) %>% 
  filter( par != 'counts') %>% 
  group_by(species, model) %>% 
  mutate( type = str_extract( par , '[a-z]+')) %>% 
  mutate( par = ifelse( type == 'par' & par == 'par', 'lambda', par)) %>%
  group_by( species, model, type ) %>% 
  arrange(par) %>% 
  mutate(par = ifelse( type == 'par'  & row_number() == max(row_number()), 'tau', par)) %>% 
  mutate( par = str_replace(par, 'par', 'alpha')) %>% 
  arrange( species, model, par)


pars[[1]]
ann_plant_mod <- function(x, form, pars) { 
  lambda <- pars[1]   
  alphas <- pars[2:(length(pars)-1)]
  tau <- pars[length(pars)]
  mm <- model.matrix( form, as.data.frame( x))
  y <- lambda*(1 + mm%*%alphas)^(tau)    
  return(y)
}
lambda <- pars[[1]][1]
alphas <- pars[[1]][2:(length(pars[[1]])-1)]
tau <- pars[[1]][length(pars[[1]])]

new_data <- matrix( c(0.5, 0.5, 0.5), nrow = 1, ncol = 3)
colnames(new_data) <- c('N1', 'N2', 'N3')
ann_plant_mod(new_data, forms[[1]], pars = pars[[1]] )
