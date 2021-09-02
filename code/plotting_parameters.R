library(tidyverse)

# ggplot themes ------------------------------------------------ # 

journal_theme <- 
  theme_bw() + 
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_text(size = 15),
          axis.text = element_text( size = 12, color = 'black'), 
          axis.title = element_text( size = 15, color = 'black'), 
          legend.text = element_text( size = 12), 
          legend.title = element_text(size = 15), 
          axis.line = element_line(color = 'black'), 
          axis.ticks = element_line(color = 'black'), 
          axis.ticks.length = unit(0.25, 'line')) 

# colors and labels ----------------------------
my_colors <- c(1,2,4,5)
my_lntps <- c(1,3,2)

species_labs <- c('Early', 'Mid', 'Late')
species_labs <- factor(species_labs, levels = c('Early', 'Mid', 'Late'), ordered = T)
