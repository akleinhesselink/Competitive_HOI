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
tiny <- .Machine$double.eps
times <- 125             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.9, 2.3) # max uptake rates mm of water per g of plant per day
K <- c(110, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedlings <- c(1, 0, 0)      # number of seedlings 
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)

parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q)

plot_transpiration(parms, my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

State <- c(soil_m, seedlings*seedling_mass)
out <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
plot_timeseries(out, parms, col = my_colors)
par(mfrow = c(1,1))
plot(out[,1:2], type = 'l', ylab = 'Soil moisture', xlab = 'day', ylim = c(0, 200))
plot(out[,c(1,3)], type = 'l', ylab = 'Biomass', xlab = 'day')

seeds <- c(1,1,1)
State <- c(soil_m, seeds*seedling_mass)
out <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )

plot(out)
plot_timeseries(out, parms, col = my_colors)
# -------- simulate annual plant experiments -------------------------- # 
comp_grad <- c(0:13)  # number of competitors
experiments <- expand.grid(as.list(rep(list ( comp_grad), 3))) # response surface experiment 
experiments <- experiments[-1, ]

experiments <- experiments %>% distinct()

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

nms <- c('focal', 'C1', 'C2')

form1 <- paste('~ -1 + focal + C1 + C2')
form2 <- paste(form1, '+ I(focal*C1) + I(focal*C2) + I(C1*C2)')

all_forms <- c(form1, form2)

gg_intra <- gg_inter <- gg_phenology <- list()
results <- est_pars <- est_alpha <- est_alpha2 <- est_alpha3 <- est_alphaHOI <- est_alphaHOI2 <- list()
i = 1

# fit annual plant model parameters --------------------------------------- # 
for(i in 1:3){ 
  # loop through each species and fit annual plant model parameters ------------------------------------------- # 
  data1 <- data[data[,i] > 0, ]  # select focal species
  data1[, i] <- data1[, i] - 1            # remove focal from competive neighborhood 
  data1$y  <- data1[, i + 3]              # focal per capita seed production
  data1$data <- as.matrix(data1[,c(1:3)])
    
  cmpttrs <- c('N1', 'N2', 'N3')
  names(cmpttrs)[i] <- nms[1]
  names(cmpttrs)[-i] <- nms[2:3]
  
  forms <- str_replace_all(all_forms, cmpttrs)
  terms <- lapply( forms, function(x) length(str_split(x, pattern = '\\+')[[1]] ))
  forms <- lapply(forms, as.formula)

  pars <- lapply( terms, function(x, ...) {c( ..., rep(0, x - 1), -1 )} , Z = lambda[i] )  
  out <- mapply(par = pars, form = forms, FUN = optim, MoreArgs = list(data = data1, fn = mod_inter, method = 'BFGS', control = list( maxit = 1e9)), SIMPLIFY = F)
  
  results[[i]] <- out 
  
  pars <- lapply( out , function(x) x$par)
  predictions <- mapply(par = pars, form = forms, FUN = mod_inter, MoreArgs = list(data = data1, predict = T))
  predictions <- data.frame(predictions)  
  names(predictions) <- paste0('type', 1:length(forms))
  data1 <- cbind(data1, predictions)

  data1 <- 
    data1 %>% 
    select ( -data) %>% 
    gather( type, val, c(y, starts_with('type'))) 
  
  data1$type <- factor(data1$type, labels = c('pred1', 'HOI', 'observed'))
  data1$type_form <- factor(data1$type, labels = c('inter + intra', '...+HOI', 'observed'))
  data_intra <- data1[ rowSums( data1[, c(1:3)[-i]] ) == 0, ]
  
  gg_intra[[i]] <- 
    ggplot( data_intra %>% filter( type != 'observed'), aes_string( x = names(data1)[i], y = names(data1)[i + 3])) + 
    geom_point( color = my_colors[i]) + 
    geom_line(aes( y = val), color = my_colors[i], linetype = 2) + 
    xlab( paste0( 'Number of intraspecific competitors')) +
    ylab( paste0( 'Per capita seed production')) +
    ggtitle( paste0('Species ', i)) + 
    my_theme + facet_wrap(~type_form)
  
  data_inter <- data1[ data1[ , i ] == 1, ] 
  
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

}

gg_intra[[1]] %+% subset( gg_intra[[1]]$data, type %in% c('pred1', 'HOI'))
gg_intra[[2]] %+% subset( gg_intra[[2]]$data, type %in% c('pred1', 'HOI'))
gg_intra[[3]] %+% subset( gg_intra[[3]]$data, type %in% c('pred1', 'HOI')) 

ggsave('figures/intra1.png', gg_intra[[1]] %+% subset( gg_intra[[1]]$data, type %in% c('pred1', 'pred2')), width = 5, height = 3.3 ) 
ggsave('figures/intra2.png', gg_intra[[2]] %+% subset( gg_intra[[2]]$data, type %in% c('pred1', 'pred2')), width = 5, height = 3.3 ) 
ggsave('figures/intra3.png', gg_intra[[3]] %+% subset( gg_intra[[3]]$data, type %in% c('pred1', 'pred2')), width = 5, height = 3.3 ) 

# Compare outcomes in 2 and 3 species communities -------------------------------- # 
gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,6))
gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,6))
gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,6))


ggsave( 'figures/inter1.png', gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,6)), width = 5.9, height = 3.3)
ggsave( 'figures/inter2.png', gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,6)), width = 5.9, height = 3.3)
ggsave( 'figures/inter3.png', gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,6)), width = 5.9, height = 3.3)


fits <- data.frame( do.call( rbind, ( lapply( results, function(x) unlist(lapply(x, function(x) x$value))))))
fits$species <- paste('species', 1:3)
fits <- fits %>% gather( model, 'MSE', X1:X2, -species)
fits$model_label <- factor(fits$model, labels = levels( gg_inter[[1]]$data$type_form)[1:2])
fits$model
fits$model_label
fits <- fits %>% group_by( species ) %>% mutate( MSE_rel = MSE/max(MSE))

ggplot( fits %>% filter( model == 'X2' ), aes( x = species, y = MSE_rel)) + 
  geom_bar(stat = 'identity')

fit_plot <- ggplot( fits, aes( x = model_label, y = MSE) ) + 
  geom_bar(stat= 'identity') + 
  facet_grid(species ~. , scales = 'free') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5 ), axis.title.x = element_blank()) 

fit_plot

ggsave('figures/fit_plot.png', fit_plot, width = 5, height = 3.7)

lapply( results, function(x) lapply(x, function(x) x$par ))
alphas

lambda <- lapply( lapply( results[[1]], function(x) x$par), head, 1)
alphas <- lapply( lapply( results[[1]], function(x) x$par), function(x) x[2:4])
tau <- lapply( lapply( results[[1]], function(x) x$par), tail, 1)

results

test <- unlist( results, recursive = F)



names (test ) <- sort( paste( c('N1', 'N2', 'N3'), rep(c('m1', 'm2'), 3), sep = '-'))
type <- rep( names(test), lapply(test, function(x) length(x$par)))
test <- data.frame( par = unlist( lapply( test, function(x) x$par) ))
test$label <- row.names( test) 

test %>% 
  separate(label, sep = '-', c('species', 'model')) %>% 
  mutate(model = str_sub( model, 1,2)) 

data.frame( type = type, )



rapply( results, function(x) x$par)

HOI <- lapply( lapply( results[[1]][4:6], function(x) x$par), function(x) tail(x, 2)[1])

quads <- lapply( lapply( results[[1]][c(2:3, 5:6)], function(x) x$par), function(x) head())

# ------------------------------------------------------------------------------ # 
gg_phenology[[1]]
gg_phenology[[2]]
gg_phenology[[3]]

# plot interaction terms 
forms
data.frame(type = c(1,2,3), forms = c('N1 + N2 + N3', 'N1 + N2 + N3 + N2*N3', 'N1 + N2 + N3' )
# Calculate with four competitors ------------------------------------------------------- # 

# parameterize model --------------------------------------------------------------------------------------------------- 
tiny <- .Machine$double.eps
times <- 125             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 3.5, 2.9, 2.3) # max uptake rates mm of water per g of plant per day
K <- c(110, 70, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedlings <- c(1, 0, 0, 0)      # number of seedlings 
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)

parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q)

plot_transpiration(parms, my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

State <- c(soil_m, seedlings*seedling_mass)
out <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
plot_timeseries(out, parms, col = my_colors)
par(mfrow = c(1,1))
plot(out[,1:2], type = 'l', ylab = 'Soil moisture', xlab = 'day', ylim = c(0, 200))
plot(out[,c(1,3)], type = 'l', ylab = 'Biomass', xlab = 'day')

seeds <- c(1,1,1, 1)
State <- c(soil_m, seeds*seedling_mass)
out <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )

plot(out)
plot_timeseries(out, parms, col = my_colors)
plot_timeseries
# -------- simulate annual plant experiments -------------------------- # 
comp_grad <- c(0:10)  # number of competitors
experiments <- expand.grid(as.list(rep(list ( comp_grad), 4))) # response surface experiment 
experiments <- experiments[-1, ]

experiments <- 
  rbind( experiments[apply( experiments, 1, function(x) sum(x == 0) == 2 ), ], 
         experiments[ apply(experiments,1, function(x) any(x == 1)) , ]  ) 


fecundity <- phenology <- data.frame( matrix( NA, nrow = nrow(experiments), ncol = 4))
per_capita_use <- use <- out <- list()

pb <- txtProgressBar(min = 1, max = nrow(experiments), style = 3)
for( i in 1:nrow(experiments)){
  setTxtProgressBar(pb , i)
  seedlings <- as.numeric(experiments[i,])
  State <- c(soil_m, seedlings*seedling_mass)
  out[[i]] <- ode(y=State, times = seq( 1, times, 0.1), func = grow, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
  
  use[[i]] <- matrix(NA, nrow = nrow(out[[i]]), ncol = 4)
  for( j in 1:4){ 
    use[[i]][, j] <- out[[i]][,2 + j]*f(out[[i]][, 2], parms$r[j], parms$K[j])
  }
  
  per_capita_use[[i]] <- sweep(use[[i]], 2, as.numeric(experiments[1, ]), '/')
  phenology[i,] <- apply( out[[i]][, c(3:6)], 2, find_phenology)
  max_biomass  <- apply( out[[i]][, c(3:6)], 2, max)
  fecundity[i,] <- (max_biomass*conversion)/seedling_mass
  
}

nrow(out[[1]])

lambda <- diag( as.matrix( fecundity [ apply(experiments, 1, sum) == 1, ]  )) # calculate lambdas 
y <- fecundity/experiments                                                    # calculate fecundity in all experiments 
data <- data.frame(experiments, y, phenology/10)
names(data ) <- c('N1', 'N2', 'N3', 'N4', 'Y1', 'Y2', 'Y3', 'Y4', paste0('PH', c(1:4)))

nms <- c('focal', 'C1', 'C2', 'C3')

form1 <- paste('~ -1 + focal + C1 + C2 + C3')
form2 <- paste(form1, '+ I(focal^2)')
form3 <- paste(form2, '+ I(C1^2) + I(C2^2) + I(C3^2)')
form4 <- paste(form1, '+ I(C1*C2) + I(C1*C3) + I(C2*C3)')
form5 <- paste(form2, '+ I(C1*C2) + I(C1*C3) + I(C2*C3)')

all_forms <- c(form1, form2, form3, form4, form5)

gg_intra <- gg_inter <- gg_phenology <- list()
est_pars <- est_alpha <- est_alpha2 <- est_alpha3 <- est_alphaHOI <- est_alphaHOI2 <- list()
i = 1

# fit annual plant model parameters --------------------------------------- # 
for(i in 1:4){ 
  # loop through each species and fit annual plant model parameters ------------------------------------------- # 
  data1 <- data[data[,i] > 1 & rowSums(data[, c(1:4)[-i]]) == 0, ]  # select focal species
  data2 <- data[data[,i] == 1 & apply( data[, c(1:4)[-i]], 1, function(x) any(x == 0 )), ] # focal with single competitor
  data3 <- rbind( data1, data2)
  data3[, i] <- data3[, i] - 1            # remove focal from competive neighborhood 
  data3$y  <- data3[, i + 3]              # focal per capita seed production
  data3$data <- as.matrix(data3[,c(1:4)])
  
  cmpttrs <- c('N1', 'N2', 'N3')
  names(cmpttrs)[i] <- nms[1]
  names(cmpttrs)[-i] <- nms[2:4]
  
  forms <- str_replace_all(all_forms, cmpttrs)
  terms <- lapply( forms, function(x) length(str_split(x, pattern = '\\+')[[1]] ))
  forms <- lapply(forms, as.formula)
  pars <- lapply( terms, function(x, ...) {c( ..., rep(1, x - 1), -1 )} , Z = lambda[i] )  
  
  data4 <- data[ data[,i] == 1, ]         # select focal species as single individual
  data4 <- rbind(data1, data4)
  data4[, i] <- data4[, i] - 1 
  data4$y <- data4[, i + 3]
  data4$data <- as.matrix(data4[, c(1:4)])
  
  data5 <- rbind(data3, data4)
  data_all <- list(data3, data3, data3, data5, data5) 
  
  out <- mapply(par = pars, form = forms, data = data_all, FUN = optim, MoreArgs = list(fn = mod_inter, control = list( maxit = 1e9 )), SIMPLIFY = F)
  
  pars <- lapply( out , function(x) x$par)
  predictions <- mapply(par = pars, form = forms, FUN = mod_inter, MoreArgs = list(data = data4, predict = T))
  predictions <- data.frame(predictions)  
  names(predictions) <- paste0('type', 1:length(forms))
  data4 <- cbind(data4, predictions)
  
  data4 <- 
    data4 %>% 
    select ( -data) %>% 
    gather( type, val, c(y, starts_with('type'))) 
  
  data4$type <- factor(data4$type, labels = c('pred1', 'pred2', 'pred3', 'HOI', 'HOI2', 'observed'))
  
  data_intra <- data4[ rowSums( data4[, c(1:4)[-i]] ) == 0, ]
  
  gg_intra[[i]] <- 
    ggplot( data_intra %>% filter( type != 'observed'), aes_string( x = names(data1)[i], y = names(data1)[i + 3])) + 
    geom_point( color = my_colors[i]) + 
    geom_line(aes( y = val), color = my_colors[i], linetype = 2) + 
    xlab( paste0( 'Number of intraspecific competitors')) +
    ylab( paste0( 'Per capita seed production')) +
    ggtitle( paste0('Species ', i)) + 
    my_theme + facet_wrap(~type)
  
  data_inter <- data4[ data4[ , i ] == 0, ] 
  
  data_inter$c1 <- data_inter[, c(1:4)[-i][1]]
  data_inter$c2 <- as.factor(data_inter[, c(1:4)[-i][2]])
  
  gg_inter[[i]] <- 
    ggplot(data_inter %>% filter( type != 'observed'), aes_string( x = 'c1', y = names(data_inter)[4:7][i], color = 'c2')) + 
    geom_point() + 
    geom_line(aes(y = val), linetype = 2) + 
    scale_color_discrete(paste( 'Density of', names(data_inter)[1:4][-i][2])) + 
    xlab( paste( 'Density of', names(data_inter)[1:4][-i][1])) + 
    ylab( paste( 'Per capita seed production of', names(data_inter)[i])) + 
    my_theme + 
    facet_wrap(~type)
  
}

gg_intra[[1]]
gg_intra[[2]]
gg_intra[[3]]

# Compare outcomes in 2 and 3 species communities -------------------------------- # 

gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,6))
gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,6))
gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,6))

# ------------------------------------------------------------------------------ # 
gg_phenology[[1]]
gg_phenology[[2]]
gg_phenology[[3]]











data4 %>% mutate( error = (Y3 - val)^2 ) %>% group_by(type) %>% summarise(MSE = mean(error, na.rm = T))




t1 <- subset(gg_inter[[1]]$data , type == 'observed' & (N3 == 0 | N2 == 0)) %>% select( val, N2, N3) 

plot_frame <- rbind( t1 %>% 
                       select(val, N2) %>% 
                       gather( competitor, density, N2) %>% 
                       filter( row_number() < 10), 
                     t1 %>% 
                       select(val, N3) %>% 
                       gather( competitor, density, N3) %>% 
                       filter( row_number() %in% c(1, 10:17)))

ggplot( plot_frame, aes(x = density, y = val, color = competitor)) + 
  geom_point() + 
  geom_line() + 
  xlab( 'Density of competitors') + 
  ylab( paste0( 'Per capita seed production of N1')) + 
  scale_color_manual(values = my_colors[c(2,3)]) + 
  my_theme


t2 <- subset(gg_inter[[2]]$data , type == 'observed' & (N1 == 0 | N3 == 0)) %>% select( val, N1, N3) 
plot_frame <- rbind( t2 %>% select(val, N1) %>% gather( competitor, density, N1) %>% filter( row_number() < 10), 
                     t2 %>% select(val, N3) %>% gather( competitor, density, N3) %>% filter( row_number() %in% c(1, 10:17)) )

ggplot( plot_frame, aes(x = density, y = val, color = competitor)) + 
  geom_point() + 
  geom_line() + 
  xlab( 'Density of competitors') + 
  ylab( paste0( 'Per capita seed production of N2')) + 
  scale_color_manual(values = my_colors[c(1,3)]) + 
  my_theme


t3 <- subset(gg_inter[[3]]$data , type == 'observed' & (N1 == 0 | N2 == 0)) %>% select( val, N1, N2) 
plot_frame <- rbind( t3 %>% select(val, N1) %>% gather( competitor, density, N1) %>% filter( row_number() < 10), 
                     t3 %>% select(val, N2) %>% gather( competitor, density, N2) %>% filter( row_number() %in% c(1, 10:17)) )

ggplot( plot_frame, aes(x = density, y = val, color = competitor)) + 
  geom_point() + 
  geom_line() + 
  xlab( 'Density of competitors') + 
  ylab( paste0( 'Per capita seed production of N3')) + 
  scale_color_manual(values = my_colors[c(1,2)]) + 
  my_theme



# re-fit with fewer 3 to 8 competitors --------------------------- 
gg_intra <- gg_inter <- gg_phenology <- list()
est_pars <- est_alpha <- list()
i = 1

for(i in 1:3){ 
  data1 <- data[data[,i] > 1 & rowSums(data[, c(1:3)[-i]]) == 0, ]  # select focal species
  data2 <- data[data[,i] == 1 & apply( data[, c(1:3)[-i]], 1, function(x) any(x == 0 )), ] # focal with single competitor
  data3 <- rbind( data1, data2)
  data3[, i] <- data3[, i] - 1            # remove focal from competive neighborhood 
  data3$y  <- data3[, i + 3]              # focal per capita seed production
  data3$data <- as.matrix(data3[,c(1:3)])
  
  temp <- data3[ apply( data3$data[, c(1:3)[-i]] , 1, function(x) any( x > 3) ), ] 
  m1 <- optim(par = c(lambda[i],1,1,1,-1), mod_inter, data = temp, control = list( maxit = 1e9 ) )
  est_alpha[[i]] <- m1$par
  
  data4 <- data[ data[,i] == 1, ]         # select focal species as single individual
  data4 <- rbind(data1, data4)
  data4[, i] <- data4[, i] - 1 
  data4$y <- data4[, i + 3]
  data4$data <- as.matrix(data4[, c(1:3)])
  
  data4$predicted <- mod_inter(pars = m1$par, data4, predict = T)
  
  data4 <- 
    data4 %>% 
    select ( -data) %>% 
    gather( type, val, c(y, predicted)) 
  
  data4$type <- factor(data4$type, labels = c('predicted', 'observed'))
  
  data_inter <- data4[ data4[ , i ] == 0, ] 
  
  data_inter$c1 <- data_inter[, c(1:3)[-i][1]]
  data_inter$c2 <- data_inter[, c(1:3)[-i][2]]
  data_inter$type <- factor(data_inter$type, levels = c('observed', 'predicted'))
  
  gg_inter[[i]] <- 
    ggplot( data_inter, aes(x = c1, y = val, color = factor(c2), linetype = type) ) + 
    geom_line() + 
    #geom_point(data = data_inter %>% filter( type == 'observed')) + 
    scale_color_discrete(paste( 'Density of', names(data_inter)[1:3][-i][2])) + 
    xlab( paste( 'Density of', names(data_inter)[1:3][-i][1])) + 
    ylab( paste( 'Per capita seed production of', names(data_inter)[i])) + 
    my_theme 
}

gg_inter[[1]] %+% 
  subset( gg_inter[[1]]$data, (N2 > 3 & N3 %in% c(0, 4, 6, 8))) + xlim(3,8) + ylim( 0, 10)

gg_inter[[2]] %+%
  subset( gg_inter[[2]]$data, (N1 > 3 & N3 %in% c(0, 4, 6, 8))) + xlim(3,8)

gg_inter[[3]] %+%
  subset( gg_inter[[3]]$data, (N1 > 3 & N2 %in% c(0, 4, 6, 8))) + xlim(3,8)

# per capita resource use 
matplot (seq(1,times,by=0.1),per_capita_use[[1]] , type = 'l')
matplot ( use[[1]], type = 'l')

per_capita_resource_use <- cbind( experiments, as.matrix( do.call( rbind, lapply( use , colSums) ) )/experiments )
plot( per_capita_resource_use [1:6, 4] , type = 'l')

experiments
per_capita_resource_use

cov(use[[1]][,1:3])

experiments[ c(1,5,6,25,26, 30,31), ] 

cov(use[[6]][, 1:3])
cov(use[[26]][, 1:3])
cov(use[[30]][, 1:3])
cov(use[[31]][, 1:3])

