rm(list = ls())

library(deSolve)

source('code/simulation_functions2.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
U <- 200                  # length of simulation in days 
R0 <- 200                 # initial resource concentration 
r <- c(4.2, 2.6, 2.0)     # max uptake rates of resource per unit mass of plant per day
#K <- c(150, 30, 0.5)      # resource half-saturation constants
a <- c(200, 60, 10)
k <- c(0.001, 0.1, 0.5)
m <- 0.09                  # tissue respiration and loss rate units of mass per unit plant mass per day 
q <- 0.07                 # resource use efficiency 
seed_mass <- c(0.005)     # seed/seedling mass 
conversion <- 0.1         # proportion live biomass converted to seed mass 

parms <- list( r = r, 
               a = a,
               k = k,
               m = m, 
               q = q, 
               R0 = R0, 
               conversion = conversion, 
               seed_mass = seed_mass, 
               #R = R, 
               U = U)

save(parms, file = 'output/parms2.rda')

# Run response surface experiments --------------------------- # 

state <- c( R0, 0, 0, 0)

B_init <- expand.grid( 
  B1 = c(0, seq(1, 8, by = 1), 16), 
  B2 = c(0, seq(1, 8, by = 1), 16), 
  B3 = c(0, seq(1, 8, by = 1), 16))

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 )  # filter out three species cases 

init <- c(10, 0, 1)
parms$n <- as.numeric(init)
state[2:4] <- as.numeric( (parms$n > 0)*parms$seed_mass )

out <- ode(state, 
                times = seq(1, parms$U, by = 0.1), 
                func = grow, 
                parms = parms, 
                rootfun = root, 
                event = list(func = event, root = T), method = 'radau')
plot(out)


out <- list() 


for( i in 1:nrow(B_init)){ 
  i <- 1
  state[2] <- parms$seed_mass
  state[3] <- parms$seed_mass
  state[4] <- parms$seed_mass 
  
  parms$n <- B_init[i, ]
  
  out[[i]] <- ode(state, 
                  times = seq(1, parms$U, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  plot(out[[1]])
}

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

save(sim_results, file = 'output/sim_results2.rda')


