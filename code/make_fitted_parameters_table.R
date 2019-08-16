rm(list = ls())
library(tidyverse)
source('code/phenomenological_models3.R')
source('code/plotting_parameters.R')

load('output/model_fits3.rda')
load('output/model_errors.rda')

model_labs <- c('Hassel', 'Hassel + HOI', 'Hassel pw', 'Model 2', 'Model 2 + HOI', 'Model 2 pw')  

pars1 <- lapply( fit1, get_fixed_pars)
pars2 <- lapply( fit2, get_fixed_pars)

pars1HOI <- lapply( fit1HOI, get_fixed_pars)
pars2HOI <- lapply( fit2HOI, get_fixed_pars)

get_pars <- function(x){ 
  data.frame( value = x, 
              par = names( x) )
}
make_df <- function(xdf){ 
  
  temp_pars <- lapply( xdf, function(x) do.call(bind_rows, lapply(x, get_pars)))
  do.call(bind_rows, 
          mapply( x = paste0('Y', 1:3), 
                  y = temp_pars, 
                  FUN = function(x, y){ 
                    y$species <- x; 
                    return(y)}, 
                  SIMPLIFY = F))
}

all_pars <- lapply( list(pars1, pars1HOI, pars2, pars2HOI), make_df)

par_table <- 
  do.call( bind_rows, 
            mapply( x = model_labs[c(1,2,4,5)], 
                    y = all_pars, 
                    FUN = function(x, y) { y$Model <- x; y }, SIMPLIFY = F)) %>% 
  separate( Model, c('Model', 'HOI'), sep = ' \\+ ', fill = 'right') %>% 
  mutate( HOI = ifelse( is.na(HOI), 'No HOI', 'HOI')) %>% 
  mutate( Model = factor( Model, levels = c('Hassel', 'Model 2'), ordered = T)) %>% 
  mutate( Species = factor( species, labels = species_labs, ordered = T)) %>% 
  mutate( HOI = factor( HOI, levels = c('No HOI', 'HOI'), ordered = T))  


par_table <- 
  par_table %>% 
  spread( par , value ) %>% 
  arrange( Model, HOI, species ) %>% 
  mutate( `Focal Species` = factor(species, labels = species_labs  )) %>% 
  left_join(error_df, by = c('Focal Species' = 'Species', 'Model', 'HOI')) %>%
  select( `Focal Species`, Model, HOI, error, lambda, contains('alpha'), tau, contains('beta') ) 


par_table %>% write_csv('output/model_parameters.csv')

