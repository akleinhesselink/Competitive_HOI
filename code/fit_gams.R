rm(list = ls())

library(mgcv)

source('code/sim_functions.R')
source('code/figure_pars.R')

data_file <- 'data/mechanistic_sim_bicultures.rds'

mydat <- readRDS(data_file)

mydat <- 
  mydat %>% 
  select(- focal_label) %>% 
  group_by( focal ) %>% 
  mutate( lambda = max(fecundity), 
          y = log(lambda/fecundity)) 

gam <- as.formula( 'y ~ te(N1) + te(N2) + te(N3)')
gam_2 <- update( gam, . ~ . + te(I(N1^2)) + te(I(N2^2)) + te(I(N3^2)))
gam_HOI <- update(gam, . ~ . + te(I(N1*N2)) + te(I(N1*N3)) + te(I(N2*N3)))

mydat <- split(mydat , mydat$focal)

fit_gam <- function( data, form, max_comp = 3){ 
  
  temp <- 
    data %>% 
    filter( comp_n < max_comp ) 
  
  gam( form, data = temp )
    
}

fits1 <- lapply(mydat, 
                FUN = fit_gam, 
                form = gam)

fits2 <- lapply(mydat, 
                FUN = fit_gam, 
                form = gam_2)

fitsHOI <- lapply(mydat, 
                FUN = fit_gam, 
                form = gam_HOI)


test <- mydat$F3
fit <- fit_gam(test, form = gam_HOI)


mydat <- mapply(fits1, mydat, FUN = function(x, y) {y$pred_1 <- predict( x, y); return(y) }, SIMPLIFY = F )
mydat <- mapply(fits2, mydat, FUN = function(x, y) {y$pred_2 <- predict( x, y); return(y) }, SIMPLIFY = F )
mydat <- mapply(fitsHOI, mydat, FUN = function(x, y) {y$pred_HOI <- predict( x, y); return(y) }, SIMPLIFY = F )

mydat <- do.call(rbind, mydat)

temp <- 
  mydat %>% 
  rename( 'obs' = y) %>% 
  gather( type, y, obs:pred_HOI) %>% 
  mutate( y = lambda/exp(y)) %>% 
  select( -lambda )

temp %>% head

temp$focal_label <- str_replace(temp$focal, 'F', 'species ')

temp %>% head

test <- temp %>% 
  spread( type, y) %>% 
  gather( type, y_pred, starts_with('pred')) %>% 
  filter( comp_n < 2) %>% 
  gather( competitor, density, starts_with('N')) %>%   
  mutate( competitor_label = str_replace(competitor, 'N', 'competitor ')) %>%
  filter( comp_n == 0 | density > 0 )

test %>% 
  ggplot(aes( x = density, y = obs, color = focal) ) + 
  geom_point(alpha = 1) + 
  facet_grid(focal_label ~ competitor_label, switch = 'both')  + 
  scale_color_manual(values = my_colors, guide = F) + 
  my_theme + 
  geom_line(aes( x = density, y = y_pred, linetype = type, group = type ))


bicultures <- 
  temp %>% 
  spread( type, y) %>% 
  gather( type, y_pred, starts_with('pred')) %>% 
  filter( comp_n  < 3 ) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  mutate( competitor_label = str_replace(competitor, 'N', 'competitor ')) %>%
  filter( comp_n == 0 | density > 0 )

bicultures %>% 
  filter( density %in% c(max(density), min(density), median(density))) %>% 
  ggplot( aes( x = density, y = obs)) + geom_point() + facet_wrap(focal_label ~ competitor_label) 

contour_data <- 
  bicultures %>% 
  mutate( z = max(obs)/obs, error = z - max(obs)/pred_1 )  %>%
  select(-competitor_label) %>% 
  spread( competitor, density, fill = 0) %>% 
  select( focal, obs, comp_n, pred_1, z, error, N1:N3) 

contour_plot <- function(cdat, focal_name){ 
  
  cdat <- cdat %>% 
    filter( focal == focal_name)
  
  p <- cdat %>% 
    filter( N3 == 0 ) %>% 
    ggplot(aes( y = N1, x = N2, z = z)) + 
    stat_contour(aes(colour = ..level..))
  
  p1 <- direct.label(p, 'bottom.pieces')
  
  p <- cdat %>% 
    filter( N2 == 0 ) %>% 
    ggplot( aes( y = N1, x = N3, z = z)) +
    stat_contour(aes(colour = ..level..))
  
  p2 <- direct.label(p, 'bottom.pieces')
  
  p <- cdat %>% 
    filter( N1 == 0 ) %>% 
    ggplot( aes( y = N2, x = N3, z = z)) +
    stat_contour(aes(colour = ..level..))
  
  p3 <- direct.label(p, 'bottom.pieces')
  
  grid.arrange(p1, p2, p3, nrow  = 1)
}

contour_plot(cdat = contour_data, focal_name = 'F1')
contour_plot(cdat = contour_data, focal_name = 'F2')
contour_plot(cdat = contour_data, focal_name = 'F3')


plot_two_sp <- function( data, focal = 'F1', C1 = 'N1', C2 = 'N2', C3 = 'N3', 
                         C2_dens = NULL ) { 
  
  data <- data[ data$focal == focal, ] 
  
  data$C1 <- as.numeric(unlist(data[, C1]))
  data$C2 <- as.numeric(unlist(data[, C2]))
  data$C3 <- as.numeric(unlist(data[, C3]))
  
  if(is.null(C2_dens)){
    C2_factor <- as.numeric( factor(data$C2)  )
    C2_dens <- floor( seq(min(C2_factor), max(C2_factor), length.out = 4))
  } 
  
  data %>%
    filter( comp_n == 0 | ( C3 == 0 )) %>%
    select(-C3) %>%
    mutate( C2 = as.factor(C2)) %>%
    filter( as.numeric(C2) %in% C2_dens ) %>%  
    ggplot( aes( y = obs, x = C1, color = C2)) +
    geom_point() +
    scale_y_continuous(name = paste0('N', str_extract(focal, '\\d'), ' fecundity'), breaks = basic_breaks()) +
    scale_x_continuous(name = paste(C1, 'density')) + 
    scale_color_discrete(name = paste(C2, 'density'))
}

test <- bicultures %>% 
  select(-competitor_label ) %>% 
  spread( competitor, density , fill = 0) %>% 
  gather( type, pred, starts_with('pred')) %>% 
  filter( type != 'pred_2')

# F1 ------------------------------------------------------------ # 
plot_two_sp(data = test, focal = 'F1', 'N1', 'N2', 'N3') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

plot_two_sp(data = test, focal = 'F1', 'N2', 'N3', 'N1') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

plot_two_sp(data = test, focal = 'F1', 'N3', 'N1', 'N2') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

# F2 ------------------------------------------------------------- # 

plot_two_sp(data = test, focal = 'F2', 'N1', 'N2', 'N3') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

plot_two_sp(data = test, focal = 'F2', 'N2', 'N3', 'N1') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

plot_two_sp(data = test, focal = 'F2', 'N3', 'N1', 'N2') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

# F3  ------------------------------------------------------------------- # 

test <- bicultures %>% 
  select(-competitor_label ) %>% 
  spread( competitor, density , fill = 0) %>% 
  gather( type, pred, starts_with('pred')) %>% 
  filter( type != 'pred_1')

plot_two_sp(data = test, focal = 'F3', 'N1', 'N2', 'N3') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

plot_two_sp(data = test, focal = 'F3', 'N2', 'N3', 'N1') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

plot_two_sp(data = test, focal = 'F3', 'N3', 'N1', 'N2') + 
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')



















bicultures %>% 
  filter( focal == 'F1') %>% 
  ggplot( aes( y = pred_1, x = obs, color = factor(comp_n))) + 
  geom_point(alpha = 0.5)

bicultures %>% 
  filter( focal == 'F2') %>% 
  ggplot( aes( y = pred_1, x = obs, color = factor(comp_n))) + 
  geom_point(alpha = 0.5)

bicultures %>% 
  filter( focal == 'F3') %>% 
  ggplot( aes( y = pred_1, x = obs, color = factor(comp_n))) + 
  geom_point(alpha = 0.5)

bicultures %>% 
  filter( focal == 'F3') %>% 
  ggplot( aes( y = pred_1_HOI, x = obs, color = factor(comp_n))) + 
  geom_point(alpha = 0.5)


