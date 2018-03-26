rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

results_file <- 'data/mechanistic_sim_monocultures.rds'

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
pulse <- 0               # amount of water supplied per day in mm  
rainy <- 10             # duration of rainy period 
r <- c(4.2, 2.9, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(98, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
epsilon <- 0.01         # rate of water evaporation and runoff mm per mm per day
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 500, length.out = 1000)
parms <- list( r = r, K = K, m =m , p = c(rep(pulse, rainy), rep(0, times - rainy)), epsilon = epsilon, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)
save(parms, file = 'data/mechanistic_parms.rda')
rm( list = names(parms) )

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# plot time series and save 

# -----------------------------------------------------
nspp <- length(parms$K)
maxdens <- 8
base <- 1.5 

experiments <- make_experiments(maxdens, base, nspp)

experiments <- make_monoculture(experiments)

experiments_with_focal <- 
  rbind( add_focal(experiments, c(1,0,0), focal_lab = 'F1'), 
       add_focal(experiments, c(0,1,0), focal_lab = 'F2'), 
       add_focal(experiments, c(0,0,1), focal_lab = 'F3')) 

results <- experiments_with_focal

seedling_mat <- split(as.matrix(results[1:3]), 1:nrow(results))
results[ , grep('^N', names(results)) ] <- do.call(rbind, 
                                                   mclapply(seedling_mat, 
                                                            FUN = run_experiment, 
                                                            parms = parms, 
                                                            mc.cores = 4))

results[ , grep('N', names(results)) ] <- (results %>% select(starts_with('N')))/(experiments_with_focal %>% select(starts_with('N')))

names( results ) <- str_replace( names(results), '^N', 'F')

results2 <- merge( experiments, results, by = c('id')) %>% 
  gather( focal2, fecundity, starts_with('F', ignore.case = F)) %>% 
  filter( focal == focal2) %>% 
  select(-focal2) %>%
  gather( competitor, density, starts_with('N')) %>% 
  group_by( id, focal ) %>% 
  mutate( comp_n = sum(density > 0)) %>% 
  spread( competitor, density, fill = 0) %>% 
  ungroup()

results2$focal_label <- paste0( 'focal\n', str_replace( results2$focal, 'F', 'N'))

results2 %>%
  filter( comp_n  < 2 ) %>% 
  gather( competitor, density, starts_with('N')) %>% 
  filter( comp_n == 0 | density > 0 ) %>%
  ggplot(aes( x = density, y = fecundity, color = competitor) ) + 
  geom_point(alpha = 1) + 
  geom_line(alpha = 0.5) + 
  scale_color_discrete(guide = F) + 
  facet_grid(focal_label ~ competitor)

saveRDS(results2, file = results_file )
