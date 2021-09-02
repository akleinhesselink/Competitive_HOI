rm(list = ls())

library(deSolve)
library(tidyverse)
source('code/simulation_functions.R')

# parameterize model --------------------------------------------------------------------------------------------------- 
U <- 200                    # Length of simulation in days 
R0 <- 400                   # Initial resource concentration 
Vmax <- 1                 # Max R uptake rates per unit root surface area  
p <- 0.5                    # Proportion root biomass per total biomass  
d <- c(0.06, 0.12, 0.36) # Root density (g/cm^3)
K <- 350                    # Resource half-saturation constants
nu <- 2/3                   # Scaling of root volume to root surface area
m <- c(0.3, 0.15, 0.053)   # Biomass loss rate per unit mass 
q <- 0.2                    # Resource use efficiency 
seed_mass <- c(0.005)       # Seed/seedling mass 
conversion <- 0.1           # Proportion live biomass converted to seed mass 
S <- 0


parms <- list( Vmax = Vmax, 
               K = K, 
               m = m, 
               q = q, 
               d = d, 
               nu = nu,
               p = p,
               R0 = R0,
               S = S,
               conversion = conversion, 
               seed_mass = seed_mass, 
               U = U)

save(parms, file = 'output/parms.rda')

# Run response surface experiments --------------------------- # 

state <- c( R = R0, B1 = 0, B2 = 0, B3 = 0)

B_scen <- expand.grid( B1 = c(0, 1), B2 = c(0,1), B3 = c(0,1))
B_scen <- B_scen[c(2, 3, 5, 8), ]

run_basic_scenarios <- function(scenarios, parms){ 
  out <- list() 
  
  for( i in 1:nrow(scenarios)){ 
    parms$n <- as.numeric( scenarios[i, 1:3] )
    parms$n[parms$n == 0 ] <- 1

    state[2] <- scenarios[i,1]*parms$seed_mass
    state[3] <- scenarios[i,2]*parms$seed_mass
    state[4] <- scenarios[i,3]*parms$seed_mass 

    out[[i]] <- ode(state, 
                times = seq(1, parms$U, by = 0.1), 
                func = grow, 
                parms = parms, 
                rootfun = root3, 
                event = list(func = event, root = T))
  }
  return(out)
}

basics <- run_basic_scenarios(B_scen, parms)

labels <- c('early', 'mid', 'late')
cols <- c('black', 'red', 'blue')

show_results <- function( out, labels, cols ){ 
  
  par(mfrow = c(length(out),1), mar = c(2,5,1,1))
  for(i in 1:(length(out)-1) ){ 
    plot(out[[i]][, 1], out[[i]][, i+2], type = 'l', ylab = labels[i], 
         col = cols[i], lty = 1)
  }
  par(mar = c(2,5,1,1), xpd = T)

  matplot(out[[length(out)]][, 1], 
          out[[length(out)]][, c((length(out) - 1):(length(out) + 1))], 
          type = 'l', 
          ylab = 'Biomass', 
          xlab = 'Days', col = cols, lty = 1)
  y_pos <- 0.8*max( apply( out[[4]][, 3:(length(out) + 1)], 2, max) )
  legend(x = 160, y = y_pos, legend = labels, col = cols, lty = 1, cex = 1.2)

}


show_results(basics, labels, cols)

par(mfrow = c(1,1), mar = c(4,4,3,3))
plot( d, m, type = 'b', col = cols, ylab = 'Loss rate (m)', 
      xlab = expression( Root ~ Density ~ (g/cm^3) ))
text(d[1], m[1], labels[1], adj = c(-0.5, 0.5))
text(d[2], m[2], labels[2], adj = c(-0.5, -1))
text(d[3], m[3], labels[3], adj = c(0.5, -1))

source('code/plotting_parameters.R')

plotting_df <- 
  basics[[4]] %>% 
  data.frame() %>% 
  gather( var, val, R:B3) %>% 
  mutate( type = ifelse(var == 'R', 'Resource', 'Biomass')) %>% 
  mutate( var_lab = factor(var, labels = c('Early', 'Mid', 'Late', 'Early')))
  
  
axis_df <- 
  plotting_df %>% 
  group_by( type) %>% 
  summarise( mx_y = max(val), mn_y = min(val)) %>%  
  ungroup()  %>% 
  mutate( y_pos = mn_y + (mx_y - mn_y)*0.5, 
          x_pos = -Inf) %>%
  mutate( labels = c('Biomass (g)', 'Resource (g/kg)'))

letters_df <- 
  plotting_df %>% 
  distinct( type) %>%  
  mutate( y_pos = Inf, 
          x_pos = -Inf, 
          label = paste0( LETTERS[1:2], ')' ))
  
mp <- 
  plotting_df %>% 
  ggplot(aes( x = time, y = val, color = var_lab )) + 
  geom_line() + 
  geom_text(data = letters_df, 
            aes( x = x_pos, y = y_pos, label = label), 
            color = 'black', 
            show.legend = F, 
            size = 5, 
            vjust = -1, 
            hjust = -1) + 
  facet_wrap(~type, ncol = 1, scales = 'free_y') + 
  scale_color_manual(values = my_colors[c(1:3, 1)], name = 'Species') + 
  xlab( 'Time (days)') + 
  scale_y_continuous(name = NULL) + 
  coord_cartesian(clip = 'off') + 
  journal_theme + 
  theme(legend.position = c(0.8, 0.8), 
        panel.spacing.y = unit(2, 'line'),
        legend.background = element_rect(fill = 'white', color = 'black'),
        legend.key.width = unit(2, 'line'), 
        plot.margin = margin( 2, 1.5, 1, 3, unit = 'line') )

mp 

mp <- 
  mp + 
  theme(strip.text = element_blank()) + 
  geom_text(data = axis_df, 
             aes( x = x_pos, 
                  y = y_pos, 
                  label = labels), color = 1, 
            angle = 90, 
            vjust = -3.5, 
            size = 5)

seed_label <- 
  plotting_df %>% 
  group_by( type, var_lab , var ) %>% 
  summarise( max_b = max(val), 
             pheno = time[ which.max(val) ]) %>%
  mutate( seeds = (max_b*parms$conversion)/parms$seed_mass ) %>% 
  filter( var != 'R') %>% 
  mutate( seeds = paste0( var_lab, '[', 't+1', ']==', round(seeds)))

seed_label$x_pos <- seed_label$pheno*c(0.6, 1, 1.15)
seed_label$y_pos <- seed_label$max_b*c(1.1, 1.1, 0.9)

mp %>% ggsave( filename = 'figures/example_dynamics_fig3.png', width = 8, height = 6 )
