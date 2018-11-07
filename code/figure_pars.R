library(ggplot2)

# graphics themes ------------------------------------------------ # 

my_theme <- theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

journal_theme <- 
  my_theme + 
  theme(  strip.background = element_blank(),
          strip.text = element_text(size = 15),
          axis.text = element_text( size = 12), 
          axis.title = element_text( size = 15), 
          legend.text = element_text( size = 12), 
          legend.title = element_text(size = 15)) 

my_colors <- c(1,2,4,5)

species_labs <- c('Early', 'Mid', 'Late')

