library(deSolve)
library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)


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
B <- c(0.01, 0.01, 0)
parms <- list(alpha = alpha, r = r)

out <- ode(B, times = 1:10, func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root)
out <- cbind( out, t(apply( out[, -1], 1, root, t = 1, parms = parms) ) )
matplot(out[, 2:4], type = 'l')

first_order <- function(x, a, t, dbdt, ... ){ 
  x + dbdt(t = 0 , B = x, ... )[[1]]*(t - a)
}

second_order <- function(x, a, t, dbdt, ... ){ 
  x + dbdt(t = 0, B = x, ... )[[1]]*(t-a) + }

points(rep(10,3),  first_order(B, a = 1, t = 10, lgrowth, parms ) )





########
par(mfrow = c(1,1))

#source('code/sim_functions.R')
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

# lgrowth <- function(t, B, parms){ 
#   with (parms, { 
#     out = r*B*(1 - alpha%*%B)    
#     return(list(out))
#   })
# }
# 
# root <- function(t, B, parms){ 
#   unlist(lgrowth(t, B, parms)) + 0.0001
# }
# 
# event <- function(t, B, parms){
#   terminate <- unlist(lgrowth(t, B, parms)) < 0 # logical vector of species to terminate
#   B[terminate] <- 0
#   return(B)
# }

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

experiments <- 
  rbind( experiments[apply( experiments, 1, function(x) sum(x == 0) == 2 ), ], 
       experiments[ apply(experiments,1, function(x) any(x == 1)) , ]  ) 

fecundity <- phenology <- data.frame( matrix( NA, nrow = nrow(experiments), ncol = 3))
per_capita_use <- use <- out <- list()

pb <- txtProgressBar(min = 1, max = nrow(experiments), style = 3)
for( i in 1:nrow(experiments)){
  setTxtProgressBar(pb , i)
  seedlings <- as.numeric(experiments[i,])
  B <- c(seedlings*seedling_mass)
  out[[i]] <- ode(y=B, times = seq( 1, times, 0.1), func = lgrowth, parms = parms, events = list(func = event, root = TRUE), rootfun = root )
  
  #phenology[i,] <- apply( out[[i]][, c(2:4)], 2, find_phenology)
  max_biomass  <- apply( out[[i]][, c(2:4)], 2, max)
  fecundity[i,] <- (max_biomass*conversion)/seedling_mass
}

lambda <- diag( as.matrix( fecundity [ apply(experiments, 1, sum) == 1, ]  )) # calculate lambdas 
y <- fecundity/experiments                                                    # calculate fecundity in all experiments 
data <- data.frame(experiments, y)
names(data ) <- c('N1', 'N2', 'N3', 'Y1', 'Y2', 'Y3')
data

gg_intra <- gg_inter <- gg_phenology <- list()
est_pars <- est_alpha <- list()
i = 1
nms <- c('focal', 'C1', 'C2')
library(stringr)
form1 <- paste('~ -1 + focal + C1 + C2')
form2 <- paste(form1, '+ I(focal^2)')
form3 <- paste(form2, '+ I(C1^2) + I(C2^2)')
form4 <- paste(form1, '+ I(C1*C2)')
form5 <- paste(form2, '+ I(C1*C2)')

all_forms <- c(form1, form2, form3, form4, form5)


for(i in 1:3){ 
  
  # loop through each species and fit annual plant model parameters ------------------------------------------- # 
  data1 <- data[data[,i] > 1 & rowSums(data[, c(1:3)[-i]]) == 0, ]  # select focal species
  data2 <- data[data[,i] == 1 & apply( data[, c(1:3)[-i]], 1, function(x) any(x == 0 )), ] # focal with single competitor
  data3 <- rbind( data1, data2)
  data3[, i] <- data3[, i] - 1            # remove focal from competive neighborhood 
  data3$y  <- data3[, i + 3]              # focal per capita seed production
  data3$data <- as.matrix(data3[,c(1:3)])
  
  cmpttrs <- c('N1', 'N2', 'N3')
  names(cmpttrs)[i] <- nms[1]
  names(cmpttrs)[-i] <- nms[2:3]
  
  forms <- str_replace_all(all_forms, cmpttrs)
  terms <- lapply( forms, function(x) length(str_split(x, pattern = '\\+')[[1]] ))
  forms <- lapply(forms, as.formula)
  pars <- lapply( terms, function(x, ...) {c( ..., rep(1, x - 1), -1 )} , Z = lambda[i] )  
  
  data4 <- data[ data[,i] == 1, ]         # select focal species as single individual
  data4 <- rbind(data1, data4)
  data4[, i] <- data4[, i] - 1 
  data4$y <- data4[, i + 3]
  data4$data <- as.matrix(data4[, c(1:3)])
  
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
  
  data_intra <- data4[ rowSums( data4[, c(1:3)[-i]] ) == 0, ]
  
  gg_intra[[i]] <- 
    ggplot( data_intra %>% filter( type != 'observed'), aes_string( x = names(data1)[i], y = names(data1)[i + 3])) + 
    geom_point( color = my_colors[i]) + 
    geom_line(aes( y = val), color = my_colors[i], linetype = 2) + 
    xlab( paste0( 'Number of intraspecific competitors')) +
    ylab( paste0( 'Per capita seed production')) +
    ggtitle( paste0('Species ', i)) + 
    my_theme + facet_wrap(~type)
  
  data_inter <- data4[ data4[ , i ] == 0, ] 
  
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
    facet_wrap(~type)
  
}

gg_intra[[1]]
gg_intra[[2]]
gg_intra[[3]]

gg_inter[[1]] %+% subset( gg_inter[[1]]$data, c2 %in% c(0,1,6))
gg_inter[[2]] %+% subset( gg_inter[[2]]$data, c2 %in% c(0,1,6))
gg_inter[[3]] %+% subset( gg_inter[[3]]$data, c2 %in% c(0,1,6))

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


















