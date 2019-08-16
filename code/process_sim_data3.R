rm(list = ls())

library(tidyverse) 

load( 'output/sim_results3.rda')
load( 'output/parms3.rda')

# convert biomass results into final seed production per plant 
sim_results <- 
  sim_results %>% 
  mutate( Y1 = parms$conversion*X2/parms$seed_mass/B1, 
          Y2 = parms$conversion*X3/parms$seed_mass/B2, 
          Y3 = parms$conversion*X4/parms$seed_mass/B3) %>% 
  select( - X1)

sim_results <- 
  sim_results %>% 
  gather( species, y, Y1, Y2, Y3)  %>% 
  mutate( species = factor(species)) %>% 
  rowwise() %>% 
  mutate( y = y + rnorm(1, 0, 0.001)) %>% # add noise to help with convergence 
  mutate( B1 = ifelse( str_extract(species, '\\d+') == 1 & B1 > 0, B1 - 1, B1)) %>%  # subtract one to find intraspecific competitor density 
  mutate( B2 = ifelse( str_extract(species, '\\d+') == 2 & B2 > 0, B2 - 1, B2)) %>% 
  mutate( B3 = ifelse( str_extract(species, '\\d+') == 3 & B3 > 0, B3 - 1, B3)) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) ) %>% 
  ungroup %>% 
  mutate( id = row_number()) %>% 
  filter( n_comp < 3)


# add HOI column to data ---------------- # 

sim_results <- 
  sim_results %>% 
  mutate( HOI = ifelse( n_comp > 1, 1, 0) ) %>% 
  rowwise() %>%
  mutate( h = max(B1*B2, B1*B3, B2*B3) )


# Filter out NA ------------------------- #  

sim_results <- 
  sim_results %>% 
  filter( !is.na(y))


cust_seq <- c(c(0:10)/10, seq(1, 9, by = 0.5))

# make continuous grid for predictions --- # 
pgrid <-
  expand.grid(species = factor( c('Y1', 'Y2', 'Y3')), B1 = cust_seq , B2 = cust_seq, B3 = cust_seq) %>% 
  data.frame() %>% 
  filter( ( (B1 > 0) + (B2 > 0) + (B3 > 0) ) < 3 ) %>% 
  ungroup() %>% 
  arrange(species, B1, B2, B3 ) %>% 
  rowwise() %>% 
  mutate( h = max(B1*B2, B1*B3, B2*B3) ) %>% 
  ungroup()  %>% 
  arrange( species, B1, B2, B3 ) %>% 
  mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) )

sim_results %>% 
  ungroup() %>% 
  filter(species == 'Y1', n_comp < 2) %>% 
  gather(comp, density, B1:B3 ) %>% 
  filter( (density > 0) | ( n_comp == 0) ) %>% 
  ggplot( aes( x = density, y = y, color = comp)) + geom_line()

sim_results %>% 
  ungroup() %>% 
  filter(species == 'Y2', n_comp < 2) %>% 
  gather(comp, density, B1:B3 ) %>% 
  filter( (density > 0) | ( n_comp == 0) ) %>% 
  ggplot( aes( x = density, y = y, color = comp)) + geom_line()

sim_results %>% 
  ungroup() %>% 
  filter(species == 'Y3', n_comp < 2) %>% 
  gather(comp, density, B1:B3 ) %>% 
  filter( (density > 0) | ( n_comp == 0) ) %>% 
  ggplot( aes( x = density, y = y, color = comp)) + geom_line()

sim_results %>% 
  ungroup() %>% 
  filter( species == 'Y3', n_comp < 3, B3 == 0 ) %>%
  select( - B3 , - starts_with('X')) %>% 
  filter( B2 %in% c(0, 2, 8), B1 < 15 ) %>% 
  ggplot( aes( y = y , x = B1, color = factor( B2 ))) + geom_point() + geom_line()

sim_results <- sim_results %>% select( -h)
pgrid <- pgrid %>% select(-h)

save(sim_results, pgrid, file = 'output/processed_results3.rda' )
