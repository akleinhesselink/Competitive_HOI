rm(list = ls())
par(mfrow = c(1,1))

source('code/sim_functions.R')
source('code/figure_pars.R')

results_file <- 'data/mechanistic_sim_bicultures.rds'

# parameterize model --------------------------------------------------------------------------------------------------- 
load(file = 'data/mechanistic_parms.rda')

plot_transpiration(parms,  my_colors)
plot_growth_rate(parms, my_colors)
plot_Rstar(parms, my_colors)

# -----------------------------------------------------
nspp <- length(parms$K)
maxdens <- 5
base <- 2 

experiments <- make_experiments(maxdens, base, nspp)

bicultures <- make_biculture(experiments)

experiments <- 
  rbind( add_focal(bicultures, c(1,0,0), focal_lab = 'F1'), 
         add_focal(bicultures, c(0,1,0), focal_lab = 'F2'), 
         add_focal(bicultures, c(0,0,1), focal_lab = 'F3')) 

results <- experiments
for( i in 1:nrow(results)){ 
  results[i, grep('N', names(results))] <- run_experiment(results[i, grep('N', names(results))], parms)
}

# rs <- split( rs, 1:nrow(rs))
# rs_results <- mclapply( rs, run_experiment, parms, mc.cores = 4)

results[ , grep('N', names(results)) ] <- (results %>% select(starts_with('N')))/(experiments %>% select(starts_with('N')))

names( results ) <- str_replace( names(results), '^N', 'F')

results2 <- merge( bicultures, results, by = c('id')) %>% 
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
