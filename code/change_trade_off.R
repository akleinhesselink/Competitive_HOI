rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
source('code/model_functions.R')
source('code/plotting_functions.R')

# new functions to adjust the parameters controlling resource uptake 
trade_off <- function(r, int = -81.8, m1 = 23.190, m2 = 7.619){ 
  new_K <- -81.8 + 23.190*r + 7.619*r^2
  return(new_K)
}

get_r <- function(step, midpoint = 2.6, upperstepsize = 0.32, lowerstepsize  = 0.1) { 
  r1 <- midpoint + step*upperstepsize    
  r3 <- midpoint - step*lowerstepsize
  r2 <- midpoint
  return ( c( r1, r2, r3 ))
}

# ---------- load original parameters ------------------- # 

load('data/parms.rda')

r <- c(2.1, 2.6, 4.2)   # max uptake rates mm of water per g of plant per day
K <- c(0.5, 30, 150)    # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  

resource_curves <- my_result <- list ()


for( k in 1:5 ) { 
  # parameterize model --------------------------------------------------------------------------------------------------- 
  r <- get_r(step = k)
  K <- trade_off(r)
  
  parms$r <- r
  parms$K <- K 
  
  resource_curves[[k]] <- plot_resource_uptake(parms, spec_labs = species_labs, R = seq(0, 200, by = 0.01))
  
  R <- 0:200

  curves <- data.frame(R = R,  mapply(x = as.list(parms$r), y = as.list(parms$K), FUN = function(x, y) { f(R = R, x, y) }) )
  
  curves <- 
    curves %>% 
    gather( species, uptake, starts_with('X')) %>% 
    mutate( species = factor(species, labels = species_labs))
  
  resource_curves[[k]] <- 
    resource_curves[[k]] + 
    journal_theme + 
    theme(legend.position = c(0.8, 0.3)) + 
    ylab('Resource uptake rate')
  
  R_init <- 200 
  seeds_init <- c(1,1,1)
  
  state <- c( R_init, NA, NA, NA)
  
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
    mutate( Competitor = factor(comp, label = species_labs)) %>% 
    mutate( species_lab = factor(species, labels = species_labs ) ) 
  
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
    journal_theme
  
  # modified Hassel with one competitor at a time --------------------------- # 
  
  # species 1 -------------------------# 
  
  form_1 <- 'y ~ lambda/(1 + (alpha.[1]*B2)^tau.[1] + (alpha.[2]*B3)^tau.[2] )'

  pw_comp_df1 <- pw_comp_df %>% 
    filter( species == 'Y1') %>% 
    spread(comp, density, fill = 0)
  
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
  
  pw_comp_df2 <- pw_comp_df %>% 
    filter( species_lab != Competitor) %>% 
    filter( species == 'Y2') %>% 
    spread(comp, density, fill = 0)
  
  
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
  
  pw_comp_df3 <- pw_comp_df %>% 
    filter( species_lab != Competitor) %>% 
    filter( species == 'Y3') %>% 
    spread(comp, density, fill = 0)
  
  
  nls3 <- 
    nls( formula = form_1, 
         data = pw_comp_df3 , 
         start = list(alpha. = c(1, 0.7), tau. = c(1, 1)), 
         algorithm = 'port', 
         lower = 0, 
         upper = 5)
  
  pw_comp_df3$pred_1 <- predict( nls3 )
  
  # try prediction on two species communities
  
  two_sp_df <- 
    results %>% 
    ungroup() %>%
    filter( n_comp < 3) %>% 
    mutate( species_lab = factor(species, labels = c('A) Early', 'B) Mid', 'C) Late')))
  
  two_sp_df$pred_y <- NA
  two_sp_df$pred_y[two_sp_df$species == 'Y1'] <- predict(nls1, newdata = two_sp_df[two_sp_df$species == 'Y1' ,] )
  two_sp_df$pred_y[two_sp_df$species == 'Y2'] <- predict(nls2, newdata = two_sp_df[two_sp_df$species == 'Y2' ,] )
  two_sp_df$pred_y[two_sp_df$species == 'Y3'] <- predict(nls3, newdata = two_sp_df[two_sp_df$species == 'Y3' ,] )
  
  my_result[[k]] <- two_sp_df
}

# plot deviations from additivity for two species competition --------------------------- # 

error_plots <- list()

for( k in 1:length(my_result)){ 

  two_sp_df <- my_result[[k]]

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
    mutate( species_lab = factor( species, labels = species_labs)) %>%
    ggplot( aes( x = species_lab, y = MSE_change, color = species_labs)) +
    geom_bar( stat = 'identity', fill = 'gray') +
    scale_color_manual(values = my_colors[1:3])  +
    ylab( 'Increase in root mean squared error') +
    xlab( 'Species') +
    journal_theme +
    theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) +
    guides( color = F)  +
    annotate( geom = 'text', 0.6, 2.7, label = 'A)', size = 5)
  
  
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
    mutate( species_lab = factor( species, labels = species_labs)) %>%
    ggplot( aes( x = species_lab, y = ME_change, color = species_labs)) +
    geom_bar( stat = 'identity', fill = 'gray') +
    scale_color_manual(values = my_colors[1:3])  +
    ylab( error_y_lab) +
    xlab( 'Species') +
    journal_theme +
    theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) +
    guides( color = F) +
    annotate( geom = 'text', 0.6, 3.65, label = 'B)', size = 5)
  
  error_plots[[k]] <- grid.arrange(MSE_plot, mean_error_plot, nrow = 1, widths = c(0.49, 0.51))
}

rcurves_data <- mapply( x = 1:length(resource_curves), y = resource_curves, function(x, y){ df <- y$data ; df$step <- x; return(df) }, SIMPLIFY = F)
rcurves_data <- do.call(rbind, rcurves_data)

my_results <- mapply( x = 1:length(resource_curves), y = my_result, function(x, y) { y$step <- x; return(y) } , SIMPLIFY = F)
my_results <- do.call(rbind, my_results)

rcurves_data <- 
  rcurves_data %>% 
  mutate( Scenario = factor( step, labels = paste( 'Scenario', 1:5))) %>% 
  rename( "Species" = species )

gg_rcurves <- 
  ggplot( data = rcurves_data, aes( x = R, y = uptake, color = Species)) + 
  geom_line() + 
  journal_theme + 
  theme( axis.text.y = element_blank()) + 
  scale_color_manual(values = my_colors) + 
  scale_y_continuous(Resource~Uptake~Rate) + 
  scale_x_continuous('Resource Availability')

gg_curves <- 
  gg_rcurves + 
  scale_x_continuous(breaks = c(0, 100, 200), 'Resource Concentration') + 
  facet_grid( ~ Scenario  )

r <- get_r(step = 1)
K <- trade_off(r)

par_tab <- data.frame( step = 1:5 ) 


MSE_plot <- 
  my_results %>%
  select( -time , -c(X2:X4), -lambda) %>%
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>%
  filter( n_comp > 0 ) %>%
  group_by(step, species, HOI) %>%
  summarise( MSE = sqrt(mean( (pred_y - y)^2)) ) %>%
  spread( HOI, MSE) %>%
  mutate( MSE_change = (`1` - `0`)  ) %>%
  ungroup() %>%
  mutate( species_lab = factor( species, labels = species_labs)) %>%
  ggplot( aes( x = step, y = MSE_change, color = species_lab)) +
  geom_point() +
  geom_line()  + 
  ylab( "Increase in root mean squared error") +
  xlab( 'Scenario') +
  journal_theme +
  theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) +
  theme(legend.position = c(0.25, 0.6)) + 
  scale_color_manual(values = my_colors[1:3])  +
  guides(color = F) + 
  ggtitle('B)') + theme(title = element_text(size = 10))


error_y_lab <- formula( Average~HOI~effect~(obs. - pred.))

mean_error_plot <- 
  my_results %>%
  select( -time , -c(X2:X4), -lambda) %>%
  filter( (B1 == 0  & species == 'Y1') | (B2 == 0 & species == 'Y2') | (B3 == 0 & species == 'Y3') )  %>%
  filter( n_comp > 0 ) %>%
  group_by(step, species, HOI) %>%
  summarise( ME = mean( (y - pred_y) ) ) %>%
  spread( HOI, ME) %>%
  mutate( ME_change = (`1` - `0`)  ) %>%
  ungroup() %>%
  mutate( species_lab = factor( species, labels = species_labs)) %>%
  ggplot( aes( x = step, y = ME_change, color = species_lab)) +
  geom_point() +
  geom_line() + 
  scale_color_manual(values = my_colors[1:3])  +
  ylab( error_y_lab) +
  xlab( 'Scenario') +
  journal_theme +
  guides( color = F) + 
  theme(axis.text.x = element_text( size = 10), axis.title.x = element_text(size = 12)) +
  ggtitle('C)') + theme(title = element_text(size = 10))


gg_curves <- 
  gg_curves +
  ggtitle('A)') + theme(title = element_text(size = 10))


error_plots <- grid.arrange(grobs = list(gg_curves, MSE_plot, mean_error_plot), 
                     layout_matrix = rbind(c(1,1), c(2,3)), heights = c(0.4, 0.6))

ggsave( error_plots, filename = 'figures/error_plots_with_trade_off.png', width = 8.5, height = 8)


# write table with parameters for Appendix ------------- # 
data.frame( do.call( cbind, par_tab %>% get_r)) %>% 
  gather( type, r, step:step.1) %>% 
  mutate( Species = factor(type, labels = species_labs )) %>% 
  mutate( Species = factor(Species, levels = species_labs, ordered = T)) %>% 
  select(Species, r) %>%
  mutate( K = trade_off(r)) %>% 
  group_by( Species) %>% 
  mutate( Scenario = row_number()) %>% 
  arrange( Scenario, Species ) %>% 
  select( Scenario, Species, r, K) %>% 
  write_csv('output/scenario_parameter_table.csv')

dev.off()
