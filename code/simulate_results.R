rm(list = ls())

library(deSolve)
library(tidyverse)
library(stringr)
library(grid)

source('code/plotting_functions.R')
source('code/model_functions.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
times <- 200             # length of simulation in days 
soil_m <- 200            # initial soil moisture (mm water in upper 500 mm of soil)
r <- c(4.2, 2.6, 2.1) # max uptake rates mm of water per g of plant per day
K <- c(150, 30, 0.5)      # resource half-saturation constant, soil moisture in mm water per 500 mm soil when plant growth is half max  
m <- 0.09                # tissue respiration and loss rate g per g per day 
q <- 0.07                 # photosynthetic water use efficiency g of carbon gain per mm of water
seedling_mass <- c(0.005) # seed/seedling mass in g 
conversion <- 0.1        # proportion live biomass converted to seed mass 
R <- seq(0, 200, length.out = 1000)
parms <- list( r = r, K = K, m =m, q = q, soil_m = soil_m, conversion = conversion, seedling_mass = seedling_mass, R = R, times = times)

save(parms,file = 'data/parms.rda')

R_init <- 200 
seeds_init <- c(1,1,1)

state <- c( R_init, seeds_init[1]*seedling_mass, seeds_init[2]*seedling_mass, seeds_init[3]*seedling_mass)

out <- ode(state, times = seq(1, 200, by = 0.01), func = grow, parms = parms, 
           rootfun = root, event = list(func = event, root = T))

out <- out[seq( 1, nrow(out), by = 50), ]

png('figures/example_timeseries.png', width = 8, height = 6, units = 'in', res = 600)


tckmark <- -0.02
par(mfrow = c(2, 2), xpd = F, 
    mai = c(0.6, 1, 0.5, 0.2), 
    mgp = c(1.5,0.2,0),  
    cex = 1.25, cex.lab = 1.2) 

plot( out[, 2], type = 'l', lwd = 3, 
      xlab = 'Time (d)', 
      ylab = 'Resource (g)', axes = F, 
      xlim = c(0, 275))
axis(side = 1, at = c(0, 100, 200, 300, 400), tck = tckmark )
axis(side = 2, at = c(0,100, 200), tck = tckmark)

box()

mtext( 'A)', adj = 0, side = 3, line = 0.2, cex = 1.4)

plot( c(1,1), type = 'n',axes = F , xlab = '', ylab = '')
par(mai = c(0.1, 0.1, 0.1, 0.2))
legend(x = 0.8, y = 1.3, 
       legend = c('Early', 'Mid', 'Late'), 
       col = my_colors[1:3], lty = 1:3, lwd = 3, 
       cex = 1.4, bty = 'n', seg.len = 5.1, y.intersp = 1.5)

par(mai = c(0.8, 1, 0.1, 0.2), 
    mgp = c(1.5,0.4,0) )

matplot( out[, c(5,4,3)], 
         type = 'l', 
         col = my_colors[c(3,2,1)], 
         lwd = 3, 
         lty = c(3,2,1), 
         xlim = c(0, 275), 
         ylab = 'Biomass (g)', 
         xlab = 'Time (d)', axes = F)
axis(side = 1, at = c(0, 100, 200, 300, 400),tck = tckmark)
axis(side = 2, at = c(0,  0.5, 1.0,  1.5), tck = tckmark)
box()

mtext( 'B)', adj = 0, side = 3, line = 0.2, cex = 1.4)

plot_R <- seq(0, 200, by = 1)

r_u <- data.frame( R = plot_R )  %>% 
  mutate( u1 = f( R, parms$r[1], parms$K[1]), 
          u2 = f( R, parms$r[2], parms$K[2]), 
          u3 = f( R, parms$r[3], parms$K[3]))

par(mai = c(0.8, 0.7, 0.1, 0.5), 
    mgp = c(1, 0.2, 0))

labelsY=parse(text=paste("Uptake~Rate~", "~'('*g~d^-1*')'", sep=""))

matplot(r_u[,1], r_u[, -1], 
        type= 'l', 
        lty = c(1:3), lwd = 3, col = my_colors,
        xlab = 'Resource (g)', 
        ylab = labelsY, 
        axes = F)

axis(side = 1, at = c(0, 100, 200), tck = tckmark)
axis(side = 2, at = c(0, 1.0, 2.0), tck = tckmark)
box()
mtext( 'C)', adj = 0, side = 3, line = 0.2, cex = 1.4)

dev.off()


# Run response surface experiments --------------------------- # 

B_init <- expand.grid( 
  B1 = c(0, seq(1, 8, by = 1), 16), 
  B2 = c(0, seq(1, 8, by = 1), 16), 
  B3 = c(0, seq(1, 8, by = 1), 16))

B_init <- B_init[-1,]

B_init <- 
  B_init %>% 
  filter( (B1 < 2) + (B2 < 2) + (B3 < 2) > 0 )  # filter out three species cases 

out <- list() 

for( i in 1:nrow(B_init)){ 
  
  state[2] <- B_init[i,1]*parms$seedling_mass
  state[3] <- B_init[i,2]*parms$seedling_mass
  state[4] <- B_init[i,3]*parms$seedling_mass 
  
  out[[i]] <- ode(state, 
                  times = seq(1, 200, by = 0.1), 
                  func = grow, 
                  parms = parms, 
                  rootfun = root, 
                  event = list(func = event, root = T), method = 'radau')
  
}

sim_results <- do.call( rbind, lapply( out, function(x) apply( x, 2, max )))
sim_results <- data.frame(B_init, sim_results )

save(sim_results, file = 'data/sim_results.rda')

dev.off()

