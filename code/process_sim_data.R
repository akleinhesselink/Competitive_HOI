
process_results <- function( x, parms, outfile = NULL){ 
  
  set.seed(0)
  
  # convert biomass results into final seed production per plant 
  sim_results <- 
    x %>% 
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
    mutate( HOI = ifelse( n_comp > 1, 1, 0) ) 
  
  # Filter out NA ------------------------- #  
  
  sim_results <- 
    sim_results %>% 
    filter( !is.na(y))
  
  max_d <- max(sim_results$B1) + 4
  
  cust_seq <- c(c(0:10)/10, seq(2, max_d, by = 0.5))
  
  # make continuous grid for predictions --- # 
  pgrid <-
    expand.grid(species = factor( c('Y1', 'Y2', 'Y3')), 
                B1 = cust_seq , B2 = cust_seq, B3 = cust_seq) %>% 
    data.frame() %>% 
    filter( ( (B1 > 0) + (B2 > 0) + (B3 > 0) ) < 3 ) %>% 
    ungroup() %>% 
    arrange(species, B1, B2, B3 ) %>% 
    ungroup()  %>% 
    arrange( species, B1, B2, B3 ) %>% 
    mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) )
  
  if(!is.null(outfile)){ 
    
    save(sim_results, pgrid, file = outfile)
    
  }else if(is.null(outfile)){ 
    
    return(sim_results)  
  }
  
  
}
  
process_results2 <- function( x, parms, outfile = NULL){ 
  
  # convert biomass results into final seed production per plant 
  # does not filter out 3-way scenarios
  sim_results <- 
    x %>% 
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
    mutate( id = row_number()) 
    #filter( n_comp < 3)
  
  
  # add HOI column to data ---------------- # 
  
  sim_results <- 
    sim_results %>% 
    mutate( HOI = ifelse( n_comp > 1, 1, 0) ) 
  
  # Filter out NA ------------------------- #  
  
  sim_results <- 
    sim_results %>% 
    filter( !is.na(y))
  
  max_d <- 40 
  cust_seq <- c(c(0:10)/10, seq(2, max_d, by = 0.5))
  #cust_seq <- sort( c( c(0:10)/10 , 2, 3, (seq(sqrt(min_d), sqrt(max_d), 1))^2) )
  
  # make continuous grid for predictions --- # 
  pgrid <-
    expand.grid(species = factor( c('Y1', 'Y2', 'Y3')), 
                B1 = cust_seq , B2 = cust_seq, B3 = cust_seq) %>% 
    data.frame() %>% 
    filter( ( (B1 > 0) + (B2 > 0) + (B3 > 0) ) < 3 ) %>% 
    ungroup() %>% 
    arrange(species, B1, B2, B3 ) %>% 
    ungroup()  %>% 
    arrange( species, B1, B2, B3 ) %>% 
    mutate( n_comp = (B1 > 0) + (B2 > 0) + (B3 > 0) )
  
  if(!is.null(outfile)){ 
    
    save(sim_results, pgrid, file = outfile)
    
  }else if(is.null(outfile)){ 
    
    return(sim_results)  
  }
}
