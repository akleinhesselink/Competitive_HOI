library(tidyverse)

x <- seq(0, 300, by = 0.1)

y <- ( 1/(1 + 1.03^(x - 100)))


resource_plot <- 
  data.frame( x = x , y = y) %>% 
  ggplot(aes( x = x, y = y)) + 
  geom_line(col = 'blue') + 
  ylab( 'Soil Resources') + 
  xlab( 'Day of Year') + 
  theme(axis.text.y = element_blank())

resource_plot

ggsave(plot = resource_plot, filename = 'ESA_2018/resource_plot1.png', width = 8, height = 5)
