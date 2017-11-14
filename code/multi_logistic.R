library(deSolve)

lgrowth <- function(t, B, parms){ 
  with (parms, { 
    out = r*B*(1 - alpha%*%B)    
    return(list(out))
  })
}

root <- function(t, B, parms){ 
  unlist(lgrowth(t, B, parms)) + 0.0001
}

event <- function(t, B, parms){
  terminate <- unlist(lgrowth(t, B, parms)) < 0 # logical vector of species to terminate
  B[terminate] <- 0
  return(B)
}


nspp <- 3 
r <- c(0.15, 0.14, 0.13)
K <- 1/r*10

alpha <- matrix(min(K), nrow = nspp, ncol = nspp)
diag(alpha) <- K
alpha <- 1/alpha
B <- c(0.01,0.01,0.01)

parms <- list(alpha = alpha, r = r)

out <- ode(B, times = 1:300, func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root)

out <- cbind( out, t(apply( out[, -1], 1, root, t = 1, parms = parms) ) )
matplot(out[, 2:4], type = 'l')


### 
r <- c(0.15, 0.14, 0.08)
K <- c(100,   200,   60)

alpha <- matrix(1e6, nrow = nspp, ncol = nspp)
diag(alpha) <- K
alpha <- 1/alpha
alpha[1,1]
alpha[1, ] <- c(alpha[1,1], 0.01, 0.01)
alpha[2, ] <- c(0.01, alpha[2,2], 0.01)
alpha[3, ] <- c(0.4e-2, 0.4e-2, alpha[3,3])

B <- rep(0.01, nspp)
parms <- list(alpha = alpha, r = r)

out <- ode(B, times = 1:300, func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root)
out <- cbind( out, t(apply( out[, -1], 1, root, t = 1, parms = parms) ) )
matplot(out[, 2:4], type = 'l')

########

library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)

par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
tiny <- .Machine$double.eps
times <- 300             # length of simulation in days 

nspp <- 3 
r <- c(0.15, 0.14, 0.13)
K <- 1/r*10
alpha <- matrix(min(K), nrow = nspp, ncol = nspp)
diag(alpha) <- K
alpha <- 1/alpha

seedlings <- c(1, 0, 0)      # number of seedlings 
seedling_mass <- c(0.01) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 

parms <- list(alpha = alpha, r = r)

B <- c(seedlings*seedling_mass)

lgrowth <- function(t, B, parms){ 
  with (parms, { 
    out = r*B*(1 - alpha%*%B)    
    return(list(out))
  })
}

root <- function(t, B, parms){ 
  unlist(lgrowth(t, B, parms)) + 0.0001
}

event <- function(t, B, parms){
  terminate <- unlist(lgrowth(t, B, parms)) < 0 # logical vector of species to terminate
  B[terminate] <- 0
  return(B)
}

out <- ode(y=B, times = seq( 1, times, 0.1), func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
matplot(out)

seeds <- c(1,1,1)
B <- c(seeds*seedling_mass)
out <- ode(y=B, times = seq( 1, times, 0.1), func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
matplot(out)

# -------- fit annual plant model parameters -------------------------- # 
comp_grad <- c(0:8)  # number of competitors
experiments <- expand.grid(as.list(rep(list ( comp_grad), 3))) # response surface experiment 
experiments <- experiments[-1, ]

fecundity <- phenology <- data.frame( matrix( NA, nrow = nrow(experiments), ncol = 3))
per_capita_use <- use <- out <- list()

pb <- txtProgressBar(min = 1, max = nrow(experiments), style = 3)
for( i in 1:nrow(experiments)){
  setTxtProgressBar(pb , i)
  seedlings <- as.numeric(experiments[i,])
  B <- c(seedlings*seedling_mass)
  out[[i]] <- ode(y=B, times = seq( 1, times, 0.1), func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
  
  phenology[i,] <- apply( out[[i]][, c(2:4)], 2, find_phenology)
  max_biomass  <- apply( out[[i]][, c(2:4)], 2, max)
  fecundity[i,] <- (max_biomass*conversion)/seedling_mass
}
phenology
find_phenology

lambda <- diag( as.matrix( fecundity [ apply(experiments, 1, sum) == 1, ]  )) # calculate lambdas 
y <- fecundity/experiments                                                    # calculate fecundity in all experiments 
data <- data.frame(experiments, y, phenology/10)
names(data ) <- c('N1', 'N2', 'N3', 'Y1', 'Y2', 'Y3', paste0('PH', c(1:3)))

gg_intra <- gg_inter <- gg_phenology <- list()
est_pars <- est_alpha <- list()
i = 1
for(i in 1:3){ 
  # loop through each species and fit annual plant model parameters ------------------------------------------- # 
  data1 <- data[data[,i] > 1 & rowSums(data[, c(1:3)[-i]]) == 0, ]  # select focal species
  data2 <- data[data[,i] == 1 & apply( data[, c(1:3)[-i]], 1, function(x) any(x == 0 )), ] # focal with single competitor
  data3 <- rbind( data1, data2)
  data3[, i] <- data3[, i] - 1            # remove focal from competive neighborhood 
  data3$y  <- data3[, i + 3]              # focal per capita seed production
  data3$data <- as.matrix(data3[,c(1:3)])
  
  m1 <- optim(par = c(lambda[i],1,1,1,-1), mod_inter, data = data3, control = list( maxit = 1e9 ) )
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
  
  data_intra <- data4[ rowSums( data4[, c(1:3)[-i]] ) == 0, ]
  
  gg_intra[[i]] <-
    ggplot( data_intra, aes_string( x = names(data1)[i], y = 'val', linetype = 'type') ) +
    geom_point(data = data_intra %>% filter( type == 'observed'), color = my_colors[i]) +
    geom_line(color = my_colors[i]) +
    scale_linetype_manual(values = c(0,2)) +
    xlab( paste0( 'Number of intraspecific competitors')) +
    ylab( paste0( 'Per capita seed production')) +
    ggtitle( paste0('Species ', i)) + 
    my_theme
  
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
  
  pheno_data <- data3
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

gg_phenology[[1]]
gg_phenology[[2]]
gg_phenology[[3]]

gg_intra[[1]]
gg_intra[[2]]
gg_intra[[3]]

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


# Compare outcomes in 2 and 3 species communities -------------------------------- # 
gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,4))
gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,4))
gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,4))


















