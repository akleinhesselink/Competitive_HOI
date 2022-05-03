

mech_model <- function( n_resource, n_species ){ 
  
  n_resource*n_species
  
}


phenom_model <- function( n_species, ord = 2 ) { 
  
  pairwise_terms <-  n_species^2 

  if( ord > 1 ) { 
    HOI_terms <- choose(n_species, ord)*( n_species - ord)
  }else{ 
    HOI_terms <- 0 
  }
  
  pairwise_terms + HOI_terms  
  
}

my_combos <- data.frame( expand.grid( n = 1:13, r = 1:13 ) )

my_combos <- my_combos[ my_combos$r == my_combos$n + 1 , ]

mech_parms <- NA
phenom_parms <- NA
phenom_parms_HOI2 <- NA

for( i in 1:nrow(my_combos)){ 
  
  mech_parms[i] <- mech_model(my_combos$r[i], my_combos$n[i])  
  phenom_parms[i] <- phenom_model(my_combos$n[i], ord = 1)
  phenom_parms_HOI2[i] <- phenom_model(my_combos$n[i], ord = 2)
}

result <- data.frame( my_combos, mech_parms, phenom_parms, phenom_parms_HOI2 )

plot(result$n, result$phenom_parms_HOI2, type = 'l', col = 'red', ylab = 'Number of parameters', xlab = 'Number of species' )
points( result$n, result$phenom_parms, type = 'l', col = 'blue')
points( result$n, result$mech_parms, type = 'l')
legend( 2, max(result$phenom_parms_HOI2)*0.9 , c('HOI model (2nd order)', 'Pairwise model', 'Mech. model (one more \nresource than species)'), lty = c(1,1,1), col = c('red', 'blue', 'black'))


