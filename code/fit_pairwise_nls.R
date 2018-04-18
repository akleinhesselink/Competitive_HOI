rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

maxiter <- 10
count <- 0 
nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'
model <- "mod_bh_ll"

basic <- 'fecundity ~ lambda.*( (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau.) '
HOI <- 'fecundity ~ lambda.*( (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 + betas.[1]*I(N1*N2) + betas.[2]*I(N1*N3) + betas.[3]*I(N2*N3))^tau.) '

lower_basic <- c(1, -1.5, rep(0, 3))
upper_basic <- c(1e3, -0.1, rep(100, 3))
lower_HOI   <- c(lower_basic, rep(-0.2, 3))
upper_HOI   <- c(upper_basic, rep(50, 3))

fits <- plot <- list()

i <- 3

mydat <- readRDS(data_file)

fit_nls <- function(data, form, init_vals, ...){
  
  temp <- 
    data %>% 
    filter(comp_n < 2)  %>%
    mutate( y = max(fecundity)/fecundity)
  
  test_fit <- nls(form, data = temp, start = init_vals, ... )
  
  return(test_fit ) 
  
} 


basic <- 'y ~ (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau. '
basic2 <- 'y ~ (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 + alpha.[4]*I(N1^2) + alpha.[5]*I(N2^2) + alpha.[6]*I(N3^2))^tau. '

type3 <- 'y ~ (1 + N1^(alpha.[1]) + N2^(alpha.[2]) + N3^(alpha.[3]))^(tau.)'

type4 <- 'y ~ 1 + (gamma.[1]*N1)^(alpha.[1]) + (gamma.[2]*N2)^(alpha.[2]) + (gamma.[3]*N3)^(alpha.[3])'

init_vals1 <- list(tau. = 1.5, alpha. = rep(1,3))
init_vals2 <- list(tau. = 1.5, alpha. = c(rep(1,3), rep(0,3)))
init_vals3 <- list(tau. = 1, alpha. = c(rep(1,3)))
init_vals4 <- list(gamma. = rep(1,3), alpha. = c(rep(1,3)))

lower1 <- c(0.1, rep(0, 3))
upper1 <- c(1.5, rep(100,3))
lower2 <- c(0.1, c(rep(0, 3), rep(-0.0001,3)))
upper2 <- c(1.5, c(rep(100, 3), rep(2,3)))
lower3 <- c(0.1, rep(0.0001, 3))
upper3 <- c(1.5, rep(5, 3))
lower4 <- c(rep(0, 3), rep(1e-5, 3))
upper4 <- c(rep(2, 3), rep(2.1, 3))

mydat <- 
  mydat %>% 
  select(- focal_label) 

mydat <- split(mydat , mydat$focal)

log(mydat[[1]]$fecundity)

fits1 <- lapply(mydat, FUN = fit_nls, form = basic, init_vals = init_vals1, lower = lower1, upper = upper1, algorithm = 'port')

fits2 <- lapply(mydat, FUN = fit_nls, form = basic2, init_vals = init_vals2, lower = lower2, upper = upper2, algorithm = 'port')

fits3 <- lapply(mydat, FUN = fit_nls, form = type3, init_vals = init_vals3, lower = lower3, upper = upper3, algorithm = 'port')

fits3

fits4 <- lapply(mydat, FUN = fit_nls, form = type4, init_vals = init_vals4, lower = lower4, upper = upper4, algorithm = 'port')

mydat[[1]]
predict(fits4[[1]])

mydat <- do.call(rbind, mydat)

mydat <- mydat %>% spread(focal, fecundity)  

mydat$F1_pred <- predict(fits1[[1]], mydat)
mydat$F2_pred <- predict(fits1[[2]], mydat)
mydat$F3_pred <- predict(fits1[[3]], mydat)

mydat$F1_pred2 <- predict(fits3[[1]], mydat)
mydat$F2_pred2 <- predict(fits3[[2]], mydat)
mydat$F3_pred2 <- predict(fits3[[3]], mydat)

temp <- 
  mydat %>% 
  gather( focal, fecundity, F1:F3_pred2) %>% 
  separate( focal, into = c('focal', 'type') , sep = '_', fill = 'right' ) %>% 
  mutate( type = ifelse(is.na(type), 'obs', type)) %>%
  spread( type, fecundity) %>% 
  group_by(focal) %>% 
  mutate( pred = max(obs)/pred) %>% 
  mutate( pred2 = max(obs)/pred2)


temp$focal_label <- str_replace(temp$focal, 'F', 'species ')

temp %>% head

density_plot <- 
  temp %>%
  filter( comp_n  < 2 ) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  mutate( competitor_label = str_replace(competitor, 'N', 'competitor ')) %>%
  filter( comp_n == 0 | density > 0 ) %>%
  ggplot(aes( x = density, y = obs, color = focal) ) + 
  geom_point(alpha = 1) + 
  #geom_line(alpha = 0.5) + 
  geom_line(aes(y=pred), linetype = 2) +
  geom_line(aes(y=pred2), linetype = 3) + 
  scale_color_manual(values = my_colors, guide = F) + 
  my_theme + 
  facet_grid(focal_label ~ competitor_label, switch = 'both') + 
  ylab( 'Per capita seed production') + 
  xlab( 'Competitor density')

density_plot

ggsave(density_plot, filename = 'figures/basic_fits.png', width = 5, height = 4)
