rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

data_file <- 'data/mechanistic_sim_bicultures.rds'

basic <- 'y ~ lambda*(1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau. '

basic_HOI <- 'y ~ lambda*(1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 +
                              beta.[1]*I(N1*N2) + beta.[2]*I(N1*N3) + beta.[3]*I(N2*N3))^tau. '

type2 <- 'y ~ lambda*(1 + N1^(alpha.[1]) + N2^(alpha.[2]) + N3^(alpha.[3]))^(tau.)'

type2_HOI <- 'y ~ lambda*(1 + N1^(alpha.[1] + beta.[1]*N2 + beta.[2]*N3) + N2^(alpha.[2] + beta.[3]*N3) + N3^(alpha.[3]))^(tau.)'

init_vals1 <- list(tau. = -1, alpha. = rep(1,3))
init_vals1_HOI <- list(tau. = -1, alpha. = rep(1,3), beta. = rep(0,3))

init_vals2 <- list(tau. = -1, alpha. = c(rep(3,3)))
init_vals2_HOI <- list(tau. = -1, alpha. = c(rep(3,3)), beta. = rep(0,3))

lower1 <- c(-1.5, rep(-1e-5, 3))
upper1 <- c(-0.1, rep(100,3))

lower1_HOI <- c(-1.5, rep(-1e-5, 3), rep(0,3))
upper1_HOI <- c(-0.1, rep(100,3), rep(1000,3))

lower2 <- c(-1.2, rep(0.0001, 3))
upper2 <- c(0.1, rep(5, 3))

mydat <- readRDS(data_file)

mydat <- 
  mydat %>% 
  select(- focal_label) %>% 
  group_by( focal ) %>% 
  mutate( lambda = max(fecundity))

mydat <- split(mydat , mydat$focal)

fits1 <- lapply(mydat, 
                FUN = fit_nls, 
                form = basic, 
                init_vals = init_vals1, 
                lower = lower1, 
                upper = upper1, 
                algorithm = 'port')

fits2 <- lapply(mydat, 
                FUN = fit_nls, 
                form = type2, 
                init_vals = init_vals2, 
                lower = lower2, 
                upper = upper2, 
                algorithm = 'port')

fits1_HOI <- lapply( mydat, 
        FUN = fit_nls, 
        form = basic_HOI, 
        init_vals = init_vals1_HOI, 
        max_comp = 3, 
        lower = lower1_HOI, 
        upper = upper1_HOI, 
        algorithm = 'port')


init_vals2_HOI <- list(tau. = -1, alpha. = c(rep(3,3)), beta. = rep(1e-3,3))
type2_HOI <- "y ~ lambda*(1 + N1^alpha.[1] + N2^alpha.[2] + N3^alpha.[3] + 
                              I(N1*N2)^beta.[1] + I(N1*N3)^beta.[2] + I(N2*N3)^beta.[3])^(tau.)"

type2_HOI2 <- "y ~ lambda*(1 + N1^alpha.[1] + N2^alpha.[2] + N3^alpha.[3] + 
                              beta.[1]*I(N1*N2) + beta.[2]*I(N1*N3) + beta.[3]*I(N2*N3))^(tau.)"

fit_nls(mydat[[3]], form = type2_HOI, 
        init_vals = init_vals2_HOI, 
        max_comp = 3,
        lower = c(lower2, c(1e-3, 1e-3, 1e-3)),
        upper = c(upper2, c(2, 1, 1)),
        algorithm = 'port')

fit_nls(mydat[[3]], form = type2_HOI2, 
        init_vals = init_vals2_HOI, 
        max_comp = 3,
        lower = c(lower2, c(1e-3, 1e-3, 1e-3)),
        upper = c(upper2, c(1, 1, 1)),
        algorithm = 'port')

hybrid_HOI <- 'y ~ lambda*((1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau. + beta.[1]*I(N1*N2) + beta.[2]*I(N1*N3) + beta.[3]*I(N2*N3))'

fit_nls(mydat[[2]], form = hybrid_HOI, 
        init_vals = init_vals1_HOI, 
        max_comp = 3, 
        lower = c(lower1, c(-1, -1, -1)), 
        upper = c(upper1, c(1, 1, 1)), 
        algorithm = 'port')


fit_nls(mydat[[3]], form = hybrid_HOI, 
        init_vals = init_vals1_HOI, 
        max_comp = 3, 
        lower = c(lower1, c(-1, -1, -1)), 
        upper = c(upper1, c(1, 1, 1)), 
        algorithm = 'port')


mydat[[3]] %>% 
  filter( comp_n < 3) %>% 
  filter( N3 == 0 )


mydat <- mapply(fits1, mydat, FUN = function(x, y) {y$pred_1 <- predict( x, y); return(y) }, SIMPLIFY = F )
mydat <- mapply(fits2, mydat, FUN = function(x, y) {y$pred_2 <- predict( x, y); return(y) }, SIMPLIFY = F )
mydat <- mapply(fits1_HOI, mydat, FUN = function(x, y) {y$pred_1_HOI <- predict( x, y); return(y) }, SIMPLIFY = F )

mydat <- do.call(rbind, mydat)

mydat <- 
  mydat %>% 
  rename( 'obs' = fecundity) %>% 
  select(-lambda)

temp <- mydat 

temp$focal_label <- str_replace(temp$focal, 'F', 'species ')

monocultures <- 
  temp %>%
  filter( comp_n  < 2 ) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  mutate( competitor_label = str_replace(competitor, 'N', 'competitor ')) %>%
  filter( comp_n == 0 | density > 0 )

monocultures %>% head

monoculture_fit_plot <- 
  monocultures %>%
  ggplot(aes( x = density, y = obs, color = focal) ) + 
  geom_point(alpha = 1) + 
  #geom_line(aes(y=pred_1), linetype = 2, alpha = 0.5) +
  #geom_line(aes(y=pred_2), linetype = 3, alpha = 0.5) + 
  geom_line(aes(y=pred_1_HOI), linetype = 4) + 
  scale_color_manual(values = my_colors, guide = F) + 
  my_theme + 
  facet_grid(focal_label ~ competitor_label, switch = 'both') + 
  ylab( 'Per capita seed production') + 
  xlab( 'Competitor density')

monoculture_fit_plot

bicultures <- 
  temp %>%
  filter( comp_n  < 3 ) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  mutate( competitor_label = str_replace(competitor, 'N', 'competitor ')) %>%
  filter( comp_n == 0 | density > 0 )


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


