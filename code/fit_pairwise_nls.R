rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

data_file <- 'data/mechanistic_sim_bicultures.rds'

basic <- 'y ~ lambda*(1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3)^tau. '
basic2 <- 'y ~ (1 + alpha.[1]*N1 + alpha.[2]*N2 + alpha.[3]*N3 + alpha.[4]*I(N1^2) + alpha.[5]*I(N2^2) + alpha.[6]*I(N3^2))^tau. '

type3 <- 'y ~ lambda*(1 + N1^(alpha.[1]) + N2^(alpha.[2]) + N3^(alpha.[3]))^(tau.)'

type4 <- 'y ~ 1 + (gamma.[1]*N1)^(alpha.[1]) + (gamma.[2]*N2)^(alpha.[2]) + (gamma.[3]*N3)^(alpha.[3])'

init_vals1 <- list(tau. = -1, alpha. = rep(1,3))
init_vals2 <- list(tau. = 1.5, alpha. = c(rep(1,3), rep(0,3)))

init_vals3 <- list(tau. = -1, alpha. = c(rep(3,3)))
init_vals4 <- list(gamma. = rep(1,3), alpha. = c(rep(1,3)))

lower1 <- c(-1.5, rep(0, 3))
upper1 <- c(-0.1, rep(100,3))
lower2 <- c(0.1, c(rep(0, 3), rep(-0.0001,3)))
upper2 <- c(1.5, c(rep(100, 3), rep(2,3)))

lower3 <- c(-1.2, rep(0.0001, 3))
upper3 <- c(0.1, rep(5, 3))

lower4 <- c(rep(0, 3), rep(1e-5, 3))
upper4 <- c(rep(2, 3), rep(2.1, 3))

mydat <- readRDS(data_file)

mydat <- 
  mydat %>% 
  select(- focal_label) %>% 
  group_by( focal ) %>% 
  mutate( lambda = max(fecundity))

mydat <- split(mydat , mydat$focal)

fits1 <- lapply(mydat, FUN = fit_nls, form = basic, init_vals = init_vals1, lower = lower1, upper = upper1, algorithm = 'port')

# fits2 <- lapply(mydat, FUN = fit_nls, form = basic2, init_vals = init_vals2, lower = lower2, upper = upper2, algorithm = 'port')
# 
fits3 <- lapply(mydat, FUN = fit_nls, form = type3, init_vals = init_vals3, lower = lower3, upper = upper3, algorithm = 'port')

# 
# fits4 <- lapply(mydat, FUN = fit_nls, form = type4, init_vals = init_vals4, lower = lower4, upper = upper4, algorithm = 'port')

mydat <- mapply(fits1, mydat, FUN = function(x, y) {y$pred_1 <- predict( x, y); return(y) }, SIMPLIFY = F )
mydat <- mapply(fits3, mydat, FUN = function(x, y) {y$pred_2 <- predict( x, y); return(y) }, SIMPLIFY = F )

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
  geom_line(aes(y=pred_1), linetype = 2) +
  geom_line(aes(y=pred_2), linetype = 3) + 
  scale_color_manual(values = my_colors, guide = F) + 
  my_theme + 
  facet_grid(focal_label ~ competitor_label, switch = 'both') + 
  ylab( 'Per capita seed production') + 
  xlab( 'Competitor density')

monoculture_fit_plot


ggsave(monoculture_fit_plot, filename = 'figures/basic_fits.png', width = 5, height = 4)



