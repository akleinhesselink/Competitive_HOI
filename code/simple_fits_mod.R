rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/sim_functions.R')
source('code/figure_pars.R')

# graphics themes ------------------------------------------------ # 

journal_theme <- 
  my_theme + 
  theme( axis.title = element_text(size = 10), 
         legend.text = element_text(size = 10), 
         legend.title = element_text(size = 12), 
         strip.text = element_text(size = 12), 
         axis.text = element_text(size = 10))

# Functions ------------------------------------------------------ #

f <- function(R, r, K){ r*R/(K + R) }              # resource (water) uptake rate. Saturates at r
dBdu <- function(u, B, R, r, K, q, m) { B*(q*f(R, r, K) - m)}  # growth as a function of biomass and resource uptake
dRdu <- function(u, B, R, r, K) { - sum(B*f(R,r, K)) } # resource (water)

grow <- function(u, State, parms, ...){
  with(parms , {
    R  <- State[1]                             # resource first
    B  <- State[2:length(State)]               # biomass for each species
    dB <- dBdu(u, B, R, r, K, q, m)
    dR <- dRdu(u, B, R, r, K)
    return( list( c(dR, dB))) } )
}

root <- function(u, State, parms) with(parms, { State[1] - m*K/(q*r-m) } )

event <- function(u, State, parms) {
  with(parms, {
    terminate <- (State[1] - m*K/(q*r-m) < 0.000001) # logical vector of species to terminate
    State[2:length(State)][ terminate ] <- 0
    return(State)
  })
}

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.6, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(150, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.1         # rate of water evaporation and runoff mm per mm per day
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 200, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)

resource_curves <- plot_resource_uptake(parms, spec_labs = c('Early', 'Mid', 'Late'), R = seq(0, 200, by = 0.01))
resource_curves

R <- 0:200
spec_labs <- c('Early', 'Mid', 'Late')
curves <- data.frame(R = R,  mapply(x = as.list(parms$r), y = as.list(parms$K), FUN = function(x, y) { f(R = R, x, y) }) )

curves <- 
  curves %>% 
  gather( species, uptake, starts_with('X')) %>% 
  mutate( species = factor(species, labels = spec_labs))

resource_curves <- resource_curves + journal_theme + theme(legend.position = c(0.8, 0.3)) + ylab('Resource uptake rate')

ggsave(filename = 'figures//resource_uptake.png', resource_curves, height = 5.5, width = 6)

R_init <- 200 
seeds_init <- c(1,1,1)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass, seeds_init[3]*seedling_mass)
test <- ode( state, times = 1:200, func = grow, parms = parms)
plot(test)

out <- ode(state, times = seq(1, 200, by = 0.01), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

ts_plots <- plot_timeseries(out, sp_labs = c('Resource', '1', '2', '3'), mytheme = journal_theme + theme(legend.position = c(0.8, 0.5), axis.title = element_text(size = 24)))

g1 <- ts_plots[[1]] + annotate(geom = 'text', 10, 220, label = 'a)', size = 5)
g2 <- ts_plots[[2]] + annotate(geom = 'text', 10, 2, label = 'b)', size = 5)

g1
g2


g3 <- resource_curves + 
  journal_theme  + 
  theme(legend.position = 'top', 
        legend.direction = 'vertical', 
        axis.text = element_blank(),
        axis.title = element_text(size = 24), 
        legend.title = element_text(size = 24), 
        legend.text = element_text(size = 20), 
        legend.key.size = unit(1.9, unit = 'lines'),
        legend.margin = margin(1.5,1,1,1, unit = 'lines')) + 
  ylab('Resource\nuptake rate') + 
  annotate( geom = 'text', 10, 2.4, label = 'c)', size = 5)

g3

ts_plot <- 
  grid.arrange( 
  g1, g2, g3,
  widths = c(1,1), 
  heights = c(1,1.2),
  layout_matrix = rbind( c(1,3), 
                         c(2,3))
  )

ggsave( ts_plot, filename = 'figures//example_timeseries.png', height = 5, width = 7)


R_init <- 200 
seeds_init <- c(0,1,0)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass, seeds_init[3]*seedling_mass)

out <- ode(state, times = seq(1, 200, by = 0.01), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

par(mfrow = c(2,1))
f_lines <- mapply( x = parms$r, y = parms$K, FUN = function(x,y) { f(out[,2], r = x, K = y)})
matplot(out[, 1], f_lines, type = 'l', col = c(1,2,4))

seeds_init <- c(1,0,0)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass, seeds_init[3]*seedling_mass)

out <- ode(state, times = seq(1, 200, by = 0.01), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

f_lines <- mapply( x = parms$r, y = parms$K, FUN = function(x,y) { f(out[,2], r = x, K = y)})
matplot(out[, 1], f_lines, type = 'l', col = c(1,2,4))


# Run response surface experiments --------------------------- # 

B_init <- expand.grid( 
  B1 = c(0, seq(1, 8, by = 1), 16), 
  B2 = c(0, seq(1, 8, by = 1), 16), 
  B3 = c(0, seq(1, 8, by = 1), 16))

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 )  # filter out three species cases 

out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  state[4] <- B_init[i,3]*parms$seedling_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, 200, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  
}

results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
results <- data.frame(B_init, results )

results

pheno <- do.call(rbind, lapply( out, function(x) apply( x, 2, which.max)))
pheno <- data.frame( B_init, pheno )

results <- 
  results %>% 
  mutate( Y1 = parms$conversion*X2/parms$seedling_mass/B1, 
          Y2 = parms$conversion*X3/parms$seedling_mass/B2, 
          Y3 = parms$conversion*X4/parms$seedling_mass/B3) %>% 
  select( - X1)

results <- 
  results %>% 
  gather( species, y, Y1, Y2, Y3)  %>% 
  mutate( y = y + rnorm(1, 0, 0.001)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>% 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( B3 = ifelse( str_extract(species, '\\d+') == 3 & B3 > 0, B3 - 1, B3)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) ) %>% 
  filter( n_comp < 3)

# add HOI column to data ---------------- # 

results <- results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) )


# Filter out NA ------------------------- #  

results <- results %>% 
  filter( !is.na(y))

# Assign lambda -------------------------- # 
results <- results %>% 
  group_by( species) %>% 
  mutate( lambda = max(y))

# model properties: 

pw_comp_df <- 
  results %>% 
  ungroup() %>%
  filter( n_comp < 2) %>% 
  gather( comp, density, B1:B3) %>% 
  filter( n_comp == 0 | density > 0) %>% 
  mutate( Competitor = factor(comp, label = c('Early', 'Mid', 'Late'))) %>% 
  #mutate( comp = ifelse(n_comp == 0, str_replace(species, 'Y', 'B'), comp )) %>% 
  mutate( species_lab = factor(species, labels = c('Early', 'Mid', 'Late') ) ) 

pw_comp_df %>% 
  filter( species_lab == Competitor)

pw_comp_gg <- 
  pw_comp_df %>% 
  filter( species_lab != Competitor) %>%
  filter( density < 15) %>% 
  ggplot( aes( x = density, y = y, color = Competitor)) + 
  geom_point() + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3]) + 
  ylab( 'Per Capita Fecundity') + 
  xlab( 'Density') + 
  my_theme


pw_comp_gg

#ggsave(pw_comp_gg + journal_theme, filename = 'ESA_2018/pairwise_comp_no_line.png', width = 10, height = 5.5)

# simple beverton holt with one exponent -------------- # 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2 + alpha.[3]*B3)^tau.'

# species 1 
pw_comp_df1 <- pw_comp_df %>% 
  filter( species == 'Y1') %>% 
  spread(comp, density, fill = 0)

nls_0_1 <- 
  nls( formula = form_0, 
       data = pw_comp_df1 , 
       start = list(alpha. = c(1, 0.7,0.3), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

nls_0_1
pw_comp_df1$pred_0 <- predict( nls_0_1 )


pw_comp_df1 <- pw_comp_df %>% 
  filter( species_lab != Competitor) %>% 
  filter( species == 'Y1') %>% 
  spread(comp, density, fill = 0)

form_0 <- 'y ~ lambda/(1 + alpha.[1]*B2 + alpha.[2]*B3)^tau.'

nls_0_1 <- 
  nls( formula = form_0, 
       data = pw_comp_df1 , 
       start = list(alpha. = c(0.7,0.3), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

nls_0_1
pw_comp_df1$pred_0 <- predict( nls_0_1 )

# species 2 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B3)^tau.'

pw_comp_df2 <- pw_comp_df %>% 
  filter( species_lab != Competitor) %>% 
  filter( species == 'Y2') %>% 
  spread(comp, density, fill = 0)

nls_0_2 <- 
  nls( formula = form_0, 
       data = pw_comp_df2 , 
       start = list(alpha. = c(1,0.3), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

pw_comp_df2$pred_0 <- predict( nls_0_2 ) 

# species 3 
form_0 <- 'y ~ lambda/(1 + alpha.[1]*B1 + alpha.[2]*B2)^tau.'

pw_comp_df3 <- pw_comp_df %>% 
  filter( species_lab != Competitor) %>% 
  filter( species == 'Y3') %>% 
  spread(comp, density, fill = 0)

nls_0_3 <- 
  nls( formula = form_0, 
       data = pw_comp_df3 , 
       start = list(alpha. = c(1, 0.7), tau. = c(1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 50)

pw_comp_df3$pred_0 <- predict( nls_0_3)

# modified Hassel with one competitor at a time

form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B2)^tau.[1] + (alpha.[2]*B3)^tau.[2] )'

# species 1 
nls1 <- 
  nls( formula = form_1, 
       data = pw_comp_df1 , 
       start = list(alpha. = c(0.7, 0.3), tau. = c( 1,1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

pw_comp_df1$pred_1 <- predict( nls1 )

# species 2 ------------------------- # 

form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B3)^tau.[2] )'

nls2 <- 
  nls( formula = form_1, 
       data = pw_comp_df2 , 
       start = list(alpha. = c(1,0.3), tau. = c(1,1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 2)

nls2
pw_comp_df2$pred_1 <- predict( nls2 )

# species 3 ------------------------- # 
form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B1)^tau.[1] + (alpha.[2]*B2)^tau.[2] )'

nls3 <- 
  nls( formula = form_1, 
       data = pw_comp_df3 , 
       start = list(alpha. = c(1, 0.7), tau. = c(1, 1)), 
       algorithm = 'port', 
       lower = 0, 
       upper = 5)

nls3
pw_comp_df3$pred_1 <- predict( nls3 )

# predict comp from pairwise model --------------------------- # 
pw_comp_df <- do.call(rbind, list( pw_comp_df1 %>% gather( comp, density, starts_with('B')), 
                     pw_comp_df2 %>% gather( comp, density, starts_with('B')), 
                     pw_comp_df3 %>% gather( comp, density, starts_with('B'))))

pw_comp_df <- 
  pw_comp_df %>% 
  gather( pred_type, pred, pred_0:pred_1)

pw_comp_df <- 
  pw_comp_df %>% 
  mutate( Model = factor( pred_type, labels = c('1', '2') )) %>%
  mutate( species_lab = factor( species_lab, labels = c('a)            Early', 'b)            Mid', 'c)            Late')))


pw_comp_pred_gg <- 
  pw_comp_df %>% 
  filter( n_comp == 0 | density > 0) %>% 
  filter( density < 15) %>% 
  ggplot( aes( x = density, y = y, color = Competitor, linetype = Model)) + 
  geom_point() + 
  geom_line(aes( y = pred ), alpha = 0.5) + 
  facet_grid(~species_lab) + 
  scale_color_manual(values = my_colors[1:3]) +
  ylab( 'Per Capita Fecundity') + 
  xlab( 'Competitor Density') + 
  my_theme + 
  journal_theme + 
  guides(color = guide_legend(order = 1)) + 
  theme( legend.position = c(0.25, 0.66), 
         strip.text = element_text(hjust = 0.1), 
         legend.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), 
         legend.title.align = c(0.5))

pw_comp_pred_gg


ggsave(pw_comp_pred_gg, filename = 'figures/pairwise_comp_with_line.png', width = 7, height = 4)

# try prediction on two species communities

two_sp_df <- results %>% 
  ungroup() %>%
  filter( n_comp < 3) %>% 
  mutate( species_lab = factor(species, labels = c('a)                      Early', 'b)                      Mid', 'c)                      Late')))


two_sp_df$pred_y <- NA
two_sp_df$pred_y[two_sp_df$species == 'Y1'] <- predict(nls1, newdata = two_sp_df[two_sp_df$species == 'Y1' ,] )
two_sp_df$pred_y[two_sp_df$species == 'Y2'] <- predict(nls2, newdata = two_sp_df[two_sp_df$species == 'Y2' ,] )
two_sp_df$pred_y[two_sp_df$species == 'Y3'] <- predict(nls3, newdata = two_sp_df[two_sp_df$species == 'Y3' ,] )

p1 <- 
  two_sp_df %>% 
  filter( B2 < 15, B3 < 15) %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y1' , B1 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B2, y = y, color = factor(B3) ) ) + 
  geom_point() + 
  scale_color_manual(values = c('gray', 'darkgray', 'black'), 'Late Species') + 
  xlab('Density Mid Species') + 
  ylab( 'Per Capita Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme( legend.title = element_text(size = 14)) + 
  theme(legend.position = c(0.7, 0.8), 
        strip.text = element_text(hjust = 0.1))

p2 <- 
  two_sp_df %>% 
  filter( B1 < 15, B3 < 15) %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y2' , B2 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B3) ) ) + 
  geom_point() + 
  scale_color_manual(values = c('gray', 'darkgray', 'black'), 'Late Species') + 
  xlab('Density Early Species') + 
  ylab( 'Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme( legend.title = element_text(size = 14)) + 
  theme(legend.position = c(0.7, 0.8), axis.title.y = element_blank(), 
        strip.text = element_text(hjust = 0.1))


p3 <- 
  two_sp_df %>% 
  filter( B1 < 15, B2 < 15) %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y3' , B3 == 0 ) %>% 
  filter( B2 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B2) ) ) + 
  geom_point() + 
  scale_color_manual(values = c('gray', 'darkgray', 'black'), 'Mid Species') + 
  xlab('Density Early Species') + 
  ylab( 'Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme( legend.title = element_text(size = 14)) + 
  theme(legend.position = c(0.7, 0.8), axis.title.y = element_blank(), 
        strip.text = element_text(hjust = 0.1))

p1
p2 
p3

n_comp2_no_line <- grid.arrange(p1, p2, p3, nrow = 1, widths = c(0.32, 0.3, 0.3))

#ggsave(n_comp2_no_line, filename = 'ESA_2018/two_sp_comp_no_line.png', width = 10, height = 5.5)

p1_line <- p1 + geom_line(aes( y = pred_y), alpha = 0.5)
p2_line <- p2 + geom_line(aes( y = pred_y), alpha = 0.5)
p3_line <- p3 + geom_line(aes( y = pred_y), alpha = 0.5)

n_comp2_pairwise_fit <- grid.arrange(p1_line, p2_line, p3_line, nrow = 1, widths = c(0.32, 0.3, 0.3))


ggsave(n_comp2_pairwise_fit, filename = 'figures/two_sp_comp_pw_line.png', width = 10, height = 5.5)

# Fit HOI ------------------------------ # 
# 
# par1 <- 
#   data.frame( species = 'Y1', coef( summary( nls1 ) )) %>% 
#   mutate( par = row.names(.)) %>% 
#   select(species, Estimate, par) %>% 
#   spread( par, Estimate)
# 
# par2 <- 
#   data.frame( species = 'Y2', coef( summary( nls2 ) )) %>% 
#   mutate( par = row.names(.)) %>% 
#   select(species, Estimate, par) %>% 
#   spread( par, Estimate)
# 
# par3 <- 
#   data.frame( species = 'Y3', coef( summary( nls3 ) )) %>% 
#   mutate( par = row.names(.)) %>% 
#   select(species, Estimate, par) %>% 
#   spread( par, Estimate)
# 
# pw_pars <- bind_rows(par1, par2, par3)
# 
# two_sp_df <- left_join(two_sp_df, pw_pars, by = 'species')
# 
# two_sp_df <- 
#   two_sp_df %>% 
#   filter(B1 < 15, B2 < 15 ,B3 < 15)
# 
# form_HOI_1 <- 'y ~ lambda/(1 + (alpha.1*B1)^tau.1 + (alpha.2*B2)^tau.2 + (alpha.3*B3)^tau.3 + beta*I(B2*B3))'
# 
# nls1_HOI <- nls( form_HOI_1, 
#                  data = two_sp_df %>% 
#                    filter( species == 'Y1', B1 == 0 ), 
#                  start = list(beta = 0), 
#                  lower = 0, 
#                  upper = 1, 
#                  algorithm = 'port')
# 
# nls1_HOI
# form_HOI_2 <- 'y ~ lambda/(1 + (alpha.1*B1)^tau.1 + (alpha.2*B2)^tau.2 + (alpha.3*B3)^tau.3 + beta*I(B1*B3))'
# 
# nls2_HOI <- nls( form_HOI_2, 
#                  data = two_sp_df %>% 
#                    filter( species == 'Y2', B2 == 0 ), 
#                  start = list(beta = 0), 
#                  lower = 0, 
#                  upper = 1, 
#                  algorithm = 'port')
# 
# 
# form_HOI_3 <- 'y ~ lambda/(1 + (alpha.1*B1)^tau.1 + (alpha.2*B2)^tau.2 + (alpha.3*B3)^tau.3 + (beta*HOI) )'
# 
# 
# nls3_HOI <- nls( form_HOI_3, 
#                  data = two_sp_df %>% 
#                    filter( species == 'Y3', B3 == 0 ), 
#                  start = list(beta = -1), 
#                  lower = c(-1), 
#                  upper = c(0), 
#                  algorithm = 'port')
# nls3_HOI
# 
# 
# # Check predictions ----------------------- # 
# two_sp_df$pred_y_HOI <- NA
# two_sp_df$pred_y_HOI[two_sp_df$species == 'Y1'] <- predict(nls1_HOI, newdata = two_sp_df[two_sp_df$species == 'Y1' ,] )
# two_sp_df$pred_y_HOI[two_sp_df$species == 'Y2'] <- predict(nls2_HOI, newdata = two_sp_df[two_sp_df$species == 'Y2' ,] )
# two_sp_df$pred_y_HOI[two_sp_df$species == 'Y3'] <- predict(nls3_HOI, newdata = two_sp_df[two_sp_df$species == 'Y3' ,] )
# 
# two_sp_df <- 
#   two_sp_df %>% 
#   gather( Fit, pred, pred_y, pred_y_HOI)
# 
# two_sp_df <- 
#   two_sp_df %>% 
#   mutate( Fit = factor( Fit, labels = c('pairwise', 'HOI')))

# two_sp_df %>% 
#   filter( (species == 'Y1' & B1 == 0) | (species == 'Y2' & B2 == 0) | (species == 'Y3' & B3 == 0) ) %>% 
#   filter( n_comp == 2) %>% 
#   group_by( species, Fit) %>% 
#   summarise( MSE = mean( ((y - pred))^2) ) %>% 
#   spread(Fit, MSE) %>% 
#   mutate( MSR_reduction = (pairwise - HOI))

# deviance(nls1)
# deviance(nls1_HOI)
# 
# deviance(nls2)
# deviance(nls2_HOI)
# 
# (deviance(nls3_HOI) - deviance(nls3))/(deviance(nls3))
# (deviance(nls2_HOI) - deviance(nls2))/(deviance(nls2))
# (deviance(nls1_HOI) - deviance(nls1))/(deviance(nls1))

MSE_lab = expression( (MSE[multi] - MSE[single])/MSE[multi])


MSE_plot <- 
  two_sp_df %>% 
  select( -time , -c(X2:X4), -lambda) %>%
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>% 
  filter( n_comp > 0 ) %>% 
  group_by(species, HOI) %>% 
  summarise( MSE = sqrt(mean( (pred_y - y)^2)) ) %>%
  spread( HOI, MSE) %>% 
  mutate( MSE_change = (`1` - `0`)  ) %>% 
  ungroup() %>% 
  mutate( species_lab = factor( species, labels = c('Early', 'Mid', 'Late'))) %>% 
  ggplot( aes( x = species_lab, y = MSE_change, color = species_lab)) + 
  geom_bar( stat = 'identity', fill = 'gray') +
  scale_color_manual(values = my_colors[1:3])  + 
  ylab( 'Increase in root mean squared error') + 
  xlab( 'Species') + 
  my_theme + 
  journal_theme + 
  theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) + 
  guides( color = F)  + 
  annotate( geom = 'text', 0.6, 2.7, label = 'a)', size = 5)



# MSE_plot <- 
#   two_sp_df %>% 
#   #filter( species == 'Y1') %>%
#   filter( Fit == 'pairwise', n_comp > 0 ) %>% 
#   select( -time , -c(X2:X4), - c(lambda:tau.3)) %>% 
#   mutate( intra = (B1 > 0  & species == 'Y1' | (B2 > 0 & species == 'Y2') | (B3 > 0 & species == 'Y3')))  %>% 
#   group_by( species, intra, HOI ) %>% 
#   summarise( MSE = mean( (pred - y)^2 ) ) %>% 
#   spread( HOI, MSE) %>% 
#   mutate( MSE_change = (`1` - `0`)  ) %>% 
#   ungroup() %>% 
#   mutate( intra_lab = factor( intra, labels = c('interspecific', 'intraspecific'))) %>% 
#   mutate( species_lab = factor( species, labels = c('Species 1', 'Species 2', 'Species 3'))) %>% 
#   ggplot( aes( x = species_lab, y = MSE_change, fill = intra_lab)) + 
#   geom_bar( stat = 'identity', position = 'dodge') + 
#   scale_fill_grey( 'HOI type' ) + 
#   scale_color_manual(values = my_colors[1:3])  + 
#   ylab( 'Increase in mean squared error') + 
#   xlab( 'Species') + 
#   my_theme + 
#   journal_theme + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text( size = 12)) + 
#   guides( color = F, fill = F)

error_y_lab <- formula( Average~HOI~effect~(obs. - pred.))

mean_error_plot <- 
  two_sp_df %>% 
  select( -time , -c(X2:X4), -lambda) %>%
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>% 
  filter( n_comp > 0 ) %>% 
  group_by(species, HOI) %>% 
  summarise( ME = mean( (y - pred_y) ) ) %>%
  spread( HOI, ME) %>% 
  mutate( ME_change = (`1` - `0`)  ) %>% 
  ungroup() %>% 
  mutate( species_lab = factor( species, labels = c('Early', 'Mid', 'Late'))) %>% 
  ggplot( aes( x = species_lab, y = ME_change, color = species_lab)) + 
  geom_bar( stat = 'identity', fill = 'gray') +
  scale_color_manual(values = my_colors[1:3])  + 
  ylab( error_y_lab) + 
  xlab( 'Species') + 
  my_theme + 
  journal_theme + 
  theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) + 
  guides( color = F) + 
  annotate( geom = 'text', 0.6, 3.65, label = 'b)', size = 5)

mean_error_plot

# mean_error_plot <- 
#   two_sp_df %>% 
#   #filter( species == 'Y1') %>%
#   filter( Fit == 'pairwise', n_comp > 0 ) %>% 
#   select( -time , -c(X2:X4), - c(lambda:tau.3)) %>% 
#   mutate( intra = (B1 > 0  & species == 'Y1' | (B2 > 0 & species == 'Y2') | (B3 > 0 & species == 'Y3')))  %>% 
#   group_by( species, intra, HOI ) %>% 
#   summarise( ME = mean( (y - pred) ) ) %>% 
#   spread( HOI, ME) %>% 
#   ungroup() %>% 
#   mutate( intra_lab = factor( intra, labels = c('interspecific', 'intraspecific'))) %>% 
#   mutate( species_lab = factor( species, labels = c('Species 1', 'Species 2', 'Species 3'))) %>% 
#   ggplot( aes( x = species_lab, y = `1`, fill = intra_lab )) + 
#   geom_bar( stat = 'identity', position = 'dodge') + 
#   scale_fill_grey( 'HOI type' ) + 
#   scale_color_manual(values = my_colors[1:3]) + 
#   ylab( error_y_lab) + 
#   xlab( 'Species') + 
#   my_theme + 
#   journal_theme + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text( size = 12)) + 
#   guides( color = F)

error_plots <- grid.arrange(MSE_plot, mean_error_plot, nrow = 1, widths = c(0.49, 0.51))

ggsave( error_plots, filename = 'figures/error_plots.png', width = 10, height = 5.5)


# Impact and sensitivity niches
# impact is the partial derivative of the resource depletion given the density of species i 
# sensitivity is the partial derivative of the per capita population growth rate to resource concentration
par(mfrow = c(1,1))

test_I <- rbind( 
  expand.grid( B1 = c(1, 1.001), B2 = 0, B3 = 0 ),
  expand.grid( B1 = 0, B2 = c(1, 1.001), B3 = 0), 
  expand.grid( B1 = 0, B2 = 0, B3 = c(1, 1.001)))

test_I
out2 <- list() 
state[1] <- parms$soil_m

for( i in 1:nrow(test_I)){ 
  
  state[2] <- test_I[i,1]*parms$seedling_mass
  state[3] <- test_I[i,2]*parms$seedling_mass
  state[4] <- test_I[i,3]*parms$seedling_mass 
  
  out2[[i]] <- ode(state, 
                  times = seq(1, 200, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  
}

plot( out2[[1]][, 2], type = 'l')
out2[[3]][1, ]
out2[[4]][1, ]
out2[[6]][1, ]
test_I

I_of_1 <- (out2[[2]][, 2] - out2[[1]][,2])/(out2[[2]][1,3] - out2[[1]][1,3] )
I_of_2 <- (out2[[4]][, 2] - out2[[3]][,2])/(out2[[4]][1,4] - out2[[3]][1,4] )
I_of_3 <- (out2[[6]][, 2] - out2[[5]][,2])/(out2[[6]][1,5] - out2[[5]][1,5] )

plot( I_of_1 , type = 'l')
plot( I_of_2, type = 'l')
plot( I_of_3, type = 'l')

sensitivity <- function( R, parms, i = 1 ) { 

  with(parms, { 
    K[i]*r[i]/(R + K[i])^2 
    })
}

parms <- parms[-which(names(parms) == 'R')]

R <- 0:200 
plot( f(R, r = parms$r[1], K = parms$K[1]), type = 'l')
points( f(R, r = parms$r[2], K = parms$K[2]), type = 'l', col = 'red')
points( f(R, r = parms$r[3], K = parms$K[3]), type = 'l', col = 'blue')

plot(x = R,  y = sensitivity(R, parms, i = 1)  , type = 'l', ylim = c(0, 0.03))
points(R,  sensitivity(R, parms, i = 2)  , type = 'l', col = 'red')
points(R,  sensitivity(R, parms, i = 3)  , type = 'l', col = 'blue')

plot(I_of_2, type = 'l')
out2[[5]][1, ]

alpha_32 <- sensitivity(out2[[5]][, 2], parms, i = 3)*I_of_2
alpha_31 <- sensitivity(out2[[5]][, 2], parms, i = 1)*I_of_1

alpha_32_w_one <- sensitivity(out2[[1]][,2], parms, i = 3)*I_of_2

plot(alpha_32, type = 'l')
plot(alpha_31, type = 'l')
which.min(alpha_32)
which.min(alpha_31)
which.min(alpha_32_w_one)

mean( alpha_32[1:1069] )
mean( alpha_31[1:921] )
mean( alpha_32_w_one[1:1021])


# fit shape of HOI function 

out[[1]]

HOI_shape_df <- 
  two_sp_df %>% 
  #filter( species == 'Y1') %>%
  filter( Fit == 'pairwise', n_comp > 0 ) %>% 
  select( -time , -c(X2:X4), - c(lambda:tau.3)) %>% 
  mutate( intra = (B1 > 0  & species == 'Y1' | (B2 > 0 & species == 'Y2') | (B3 > 0 & species == 'Y3'))) %>% 
  mutate( HOI = y  - pred) 

HOI_shape_df %>% 
  filter( species == 'Y1') %>% 
  filter( intra == F) %>% 
  ggplot(aes(x = B2, y = B3, z = HOI, fill = HOI )) + 
  geom_contour() + 
  geom_raster(interpolate = T) + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')

HOI_shape_df %>% 
  filter( species == 'Y2') %>% 
  filter( intra == F) %>% 
  filter(B1 > 0 , B3 > 0) %>%
  ggplot(aes(x = B1, y = B3, fill = HOI )) + 
  geom_raster() + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')

HOI_shape_df %>% 
  filter( species == 'Y3') %>% 
  filter( intra == F) %>% 
  filter(B1 > 0 , B2 > 0) %>%
  ggplot(aes(x = B1, y = B2, fill = HOI )) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')

HOI_shape_df %>% 
  filter( species == 'Y3') %>% 
  filter( intra == F) %>% 
  mutate( B1B2 = B1*B2) %>% 
  ggplot(aes(x = B1B2, y = HOI)) + 
  geom_point()

HOI_shape_df %>% 
  filter( species == 'Y2') %>% 
  filter( intra == F) %>% 
  filter( B2 == 0) %>% 
  mutate( B3B1 = B3*B1) %>% 
  ggplot(aes(x = B3B1, y = HOI)) + 
  geom_point()




p1 <- 
  two_sp_df %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y1' , B1 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B2, y = y, color = factor(B3) , linetype = Fit) ) + 
  geom_point() + 
  scale_color_manual(values = c('gray', 'darkgray', 'black'), 'Density Spp. 3') + 
  xlab('Density Species 2') + 
  ylab( 'Per Capita Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme(legend.title = element_text(size = 14 )) + 
  theme(legend.position = c(0.7, 0.8))

p2 <- 
  two_sp_df %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y2' , B2 == 0 ) %>% 
  filter( B3 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B3), linetype = Fit) ) + 
  geom_point() + 
  scale_color_manual(values = c('gray', 'darkgray', 'black'), 'Density Spp. 3') + 
  xlab('Density Species 1') + 
  ylab( 'Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme(legend.title = element_text(size = 14 )) + 
  theme(legend.position = c(0.7, 0.8), axis.title.y = element_blank())


p3 <- 
  two_sp_df %>% 
  mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
  filter( species == 'Y3' , B3 == 0 ) %>% 
  filter( B2 %in% c(0, 2, 8)) %>% 
  ggplot( aes( x = B1, y = y, color = factor(B2), linetype = Fit ) ) + 
  geom_point() + 
  scale_color_manual(values = c('gray', 'darkgray', 'black'), 'Density Spp. 2') + 
  xlab('Density Species 1') + 
  ylab( 'Fecundity') + 
  facet_wrap(~ species_lab) + 
  my_theme + 
  journal_theme + 
  theme(legend.title = element_text(size = 14 )) + 
  guides( color = guide_legend(order = 1)) +
  theme(legend.position = c(0.7, 0.7), axis.title.y = element_blank())

p1 
p2
p3


p1 <- p1 + geom_line(aes( y= pred)) + scale_linetype(guide = F)
p2 <- p2 + geom_line(aes( y = pred)) + scale_linetype(guide = F)
p3 <- p3 + geom_line(aes( y = pred)) 

p3

p3_test <- 
          two_sp_df %>% 
          mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
          filter( species == 'Y3' , B3 == 0 ) %>% 
          ggplot( aes( x = B1, y = y, color = factor(B2), linetype = Fit ) ) + 
          geom_point() + 
          xlab('Density Species 1') + 
          ylab( 'Fecundity') + 
          facet_wrap(~ species_lab) + 
          my_theme + 
          journal_theme + 
          theme(legend.title = element_text(size = 14 )) + 
          guides( color = guide_legend(order = 1)) +
          theme(legend.position = c(0.7, 0.7), axis.title.y = element_blank())
        
        
        p3_test + geom_line(aes( y = pred)) + scale_linetype()
        
        two_sp_df %>% 
          mutate( lambda_plot  = ifelse (y == lambda, T, F)) %>% 
          filter( species == 'Y3' , B3 == 0 ) %>%
          select( y, Fit, pred) %>% 
          mutate( error = (y - pred)^2) %>% 
          group_by( Fit) %>% 
          summarise( mean(error))
        
        
        n_comp2_HOI_fit <- grid.arrange(p1, p2, p3, nrow = 1, widths = c(0.32, 0.3, 0.3))
        
        #ggsave(n_comp2_HOI_fit, filename = 'ESA_2018/two_sp_comp_HOI_line.png', width = 10, height = 5.5)
        
        